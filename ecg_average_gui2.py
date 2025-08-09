import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirnotch, filtfilt, find_peaks
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
from matplotlib.widgets import SpanSelector
from PIL import ImageGrab
import io
import win32clipboard
import win32con
from datetime import datetime
import threading
import queue
import time
from concurrent.futures import ThreadPoolExecutor

# === Параметры сигнала ===
FS = 2000.0  # Полная частота дискретизации для математики
FS_DISPLAY = 2000.0  # Полная частота для отрисовки
BYTES_ECG = 3
BYTES_AUX = 2 * 4
N_CHANNELS = 8
SAMPLE_SIZE = BYTES_ECG * N_CHANNELS + BYTES_AUX
START_DELAY = 0.5
READ_COUNT = int(100 * FS)  # Оставляем для быстрого предпросмотра

# Параметры для обработки больших файлов
CHUNK_SIZE_MS = 30000  # 30 секунд на чанк для экономии памяти
OVERLAP_MS = 2000      # 2 секунды перекрытия между чанками

# === Функции обработки ===
def open_ecg_memmap(path):
    mm = np.memmap(path, dtype=np.uint8, mode='r')
    total = mm.size // SAMPLE_SIZE
    return mm, total

def extract_channel(mm, total, ch, start, count):
    end = min(start + count, total)
    arr = np.zeros(end - start, dtype=np.int32)
    offset = start * SAMPLE_SIZE + ch * BYTES_ECG
    for i in range(end - start):
        b0, b1, b2 = mm[offset:offset+3]
        raw = (int(b0)<<16) | (int(b1)<<8) | int(b2)
        arr[i] = raw - (1<<23)
        offset += SAMPLE_SIZE
    return arr

def detect_r_peaks_streaming(mm, total, ch, fs, chunk_size=10000, overlap=2000, progress_callback=None, **kwargs):
    """
    Последовательная детекция R-пиков для больших файлов
    """
    all_peaks = []
    chunk_samples = int(chunk_size * fs / 1000)  # размер чанка в сэмплах
    overlap_samples = int(overlap * fs / 1000)   # перекрытие в сэмплах
    
    for start in range(0, total, chunk_samples - overlap_samples):
        end = min(start + chunk_samples, total)
        
        # Извлекаем чанк данных
        chunk = extract_channel(mm, total, ch, start, end - start)
        chunk = chunk.astype(np.float64) * ((2*2.4)/(2**24)) * 1e3
        
        # Применяем фильтры
        chunk = highpass_filter(chunk, fs, cutoff=0.5)
        if kwargs.get('notch_enabled', True):
            chunk = notch_filter(chunk, fs)
        
        # Детектируем пики в чанке
        peaks, _, _, _, _ = detect_r_peaks(chunk, fs, **kwargs)
        
        # Корректируем позиции пиков относительно начала файла
        peaks += start
        
        # Фильтруем пики, которые выходят за границы файла
        peaks = peaks[(peaks >= 0) & (peaks < total)]
        
        # Добавляем пики (исключаем дубликаты в области перекрытия)
        if all_peaks:
            # Исключаем пики в области перекрытия
            overlap_start = start + overlap_samples
            peaks = peaks[peaks >= overlap_start]
        
        all_peaks.extend(peaks)
        
        # Обновляем прогресс через callback
        if progress_callback:
            progress = (end / total) * 100
            progress_callback(progress)
    
    return np.array(all_peaks)

def notch_filter(data, fs, freq=50.0, Q=30.0):
    b, a = iirnotch(freq, Q, fs)
    return filtfilt(b, a, data)

def highpass_filter(data, fs, cutoff=0.5):
    """Фильтр высоких частот для устранения плавающей изолинии"""
    from scipy.signal import butter, filtfilt
    nyquist = fs / 2
    normal_cutoff = cutoff / nyquist
    b, a = butter(2, normal_cutoff, btype='high', analog=False)
    return filtfilt(b, a, data)

def downsample_for_display(data, fs_original, fs_target):
    """Понижение частоты дискретизации для отображения"""
    if fs_original <= fs_target:
        return data
    
    # Вычисляем коэффициент понижения
    downsample_factor = int(fs_original / fs_target)
    
    # Простое понижение частоты дискретизации (каждый N-й сэмпл)
    return data[::downsample_factor]

def reduce_amplitude_resolution(data, factor=5):
    """Снижение разрешения по амплитуде для ускорения отрисовки"""
    # Округляем значения до factor знаков после запятой
    return np.round(data / factor) * factor

def detect_r_peaks(ecg, fs, sensitivity=1.0, min_rr_ms=200, window_ms=120, search_window_ms=50, notch_enabled=True):
    """
    Улучшенный алгоритм Томсона (Pan-Tompkins) для детекции R-пиков в ЭКГ
    
    Параметры:
    - sensitivity: чувствительность (0.1-2.0, по умолчанию 1.0)
    - min_rr_ms: минимальное расстояние между пиками в мс (по умолчанию 200)
    - window_ms: размер окна интегрирования в мс (по умолчанию 120)
    - search_window_ms: окно поиска точного пика в мс (по умолчанию 50)
    - notch_enabled: применять ли notch фильтр (по умолчанию True)
    """
    from scipy.signal import butter, filtfilt
    
    # 1. Предварительная фильтрация - устранение плавающей изолинии
    ecg_filtered = highpass_filter(ecg, fs, cutoff=0.5)
    
    # 2. Полосовой фильтр (5-15 Гц для QRS)
    nyquist = fs / 2
    low = 5 / nyquist
    high = 15 / nyquist
    b, a = butter(2, [low, high], btype='band')
    filtered = filtfilt(b, a, ecg_filtered)
    
    # 3. Дифференцирование
    diff_signal = np.diff(filtered)
    
    # 4. Возведение в квадрат
    squared = diff_signal ** 2
    
    # 5. Интегрирование с окном
    window_size = int(window_ms * fs / 1000)
    integrated = np.convolve(squared, np.ones(window_size), mode='same')
    
    # 6. Улучшенный адаптивный порог
    # Используем медиану вместо максимума для более стабильного порога
    threshold_i1 = (0.6 / sensitivity) * np.percentile(integrated, 95)
    threshold_i2 = 0.6 * threshold_i1
    
    # Списки для хранения пиков
    r_peaks = []
    noise_peaks = []
    signal_peaks = []
    
    # Минимальное расстояние между пиками
    min_rr = int(min_rr_ms * fs / 1000)
    
    # 7. Поиск пиков с улучшенной логикой
    for i in range(2, len(integrated) - 2):
        # Проверяем локальный максимум с более строгими условиями
        if (integrated[i] > integrated[i-1] and integrated[i] > integrated[i+1] and
            integrated[i] > integrated[i-2] and integrated[i] > integrated[i+2]):
            
            if integrated[i] > threshold_i1:
                # Проверяем минимальное расстояние
                if not r_peaks or (i - r_peaks[-1]) > min_rr:
                    r_peaks.append(i)
                    signal_peaks.append(i)
            elif integrated[i] > threshold_i2:
                noise_peaks.append(i)
    
    # 8. Адаптивное обновление порогов
    if len(signal_peaks) > 8:
        signal_level = np.mean([integrated[peak] for peak in signal_peaks[-8:]])
        noise_level = np.mean([integrated[peak] for peak in noise_peaks[-8:]]) if noise_peaks else signal_level * 0.5
        
        threshold_i1 = (noise_level + 0.25 * (signal_level - noise_level)) / sensitivity
        threshold_i2 = 0.5 * threshold_i1
    
    # 9. Поиск точных позиций R-пиков в оригинальном сигнале
    final_peaks = []
    search_window = int(search_window_ms * fs / 1000)
    
    for peak in r_peaks:
        start = max(0, peak - search_window)
        end = min(len(ecg), peak + search_window)
        
        # Ищем максимум в окне
        local_max_idx = start + np.argmax(ecg[start:end])
        final_peaks.append(local_max_idx)
    
    # Возвращаем также промежуточные сигналы для визуализации
    return np.array(final_peaks), filtered, diff_signal, squared, integrated

def average_complex(ecg, peaks, fs, window_ms=700):
    half = int(window_ms * fs / 1000 / 2)
    segs = []
    for r in peaks:
        if r-half >= 0 and r+half < len(ecg):
            segs.append(ecg[r-half:r+half])
    if not segs:
        return None, None
    A = np.vstack(segs)
    mean = A.mean(axis=0)
    t = np.linspace(-window_ms/2, window_ms/2, mean.size)
    return t, mean

def average_complex_from_peaks(mm, total, peaks, fs, window_ms=700):
    """Усреднение комплекса из пиков, извлекая данные из memmap"""
    if len(peaks) < 2:
        return None, None
    
    half = int(window_ms * fs / 1000 / 2)
    complexes = []
    
    for peak in peaks:
        # Проверяем, что пик находится в допустимых границах
        if peak < 0 or peak >= total:
            continue
            
        start = max(0, peak - half)
        end = min(total, peak + half)
        
        # Проверяем, что у нас достаточно данных для комплекса
        if end - start >= half:  # Минимум половина окна
            try:
                # Извлекаем данные из memmap
                raw = extract_channel(mm, total, 0, start, end - start).astype(np.float64)
                ecg = raw * ((2*2.4)/(2**24)) * 1e3
                
                # Применяем фильтры
                ecg = highpass_filter(ecg, fs, cutoff=0.5)
                ecg = notch_filter(ecg, fs)
                
                # Дополняем до полного размера если нужно
                if len(ecg) < 2 * half:
                    # Дополняем нулями до полного размера
                    padding = 2 * half - len(ecg)
                    if start == 0:  # Дополняем справа
                        ecg = np.pad(ecg, (0, padding), mode='constant')
                    else:  # Дополняем слева
                        ecg = np.pad(ecg, (padding, 0), mode='constant')
                
                complexes.append(ecg)
            except Exception as e:
                print(f"Ошибка при обработке пика {peak}: {e}")
                continue
    
    if not complexes:
        return None, None
    
    # Приводим все комплексы к одинаковому размеру
    min_length = min(len(c) for c in complexes)
    complexes = [c[:min_length] for c in complexes]
    
    A = np.vstack(complexes)
    mean = A.mean(axis=0)
    t = np.linspace(-window_ms/2, window_ms/2, mean.size)
    return t, mean

def create_ecg_markers(t, signal):
    """
    Создание пустого словаря маркеров ЭКГ - пользователь добавит их вручную
    """
    markers = {}
    return markers

def copy_image_to_clipboard(img):
    output = io.BytesIO()
    img.convert("RGB").save(output, "BMP")
    data = output.getvalue()[14:]
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32con.CF_DIB, data)
    win32clipboard.CloseClipboard()



def on_close(app):
    # Очищаем ресурсы
    if hasattr(app, 'executor'):
        app.executor.shutdown(wait=False)
    plt.close('all')
    app.destroy()
    sys.exit(0)

# === Главное окно ===
class ECGApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ECG Analysis Interface")
        self.geometry("1200x800")
        self.protocol("WM_DELETE_WINDOW", lambda: on_close(self))
        
        # --- Меню ---
        menubar = tk.Menu(self)
        
        # Меню Файл
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Открыть запись...", command=self.open_file)
        filemenu.add_separator()
        filemenu.add_command(label="Сохранить параметры...", command=self.save_params)
        filemenu.add_command(label="Загрузить параметры...", command=self.load_params)
        filemenu.add_separator()
        filemenu.add_command(label="Выход", command=lambda: on_close(self))
        menubar.add_cascade(label="Файл", menu=filemenu)
        
        # Меню Настройки
        settingsmenu = tk.Menu(menubar, tearoff=0)
        settingsmenu.add_command(label="Параметры алгоритма Томсона...", command=self.show_tompson_settings)
        settingsmenu.add_separator()
        settingsmenu.add_command(label="Сброс параметров", command=self.reset_tompson_params)
        menubar.add_cascade(label="Настройки", menu=settingsmenu)
        
        self.config(menu=menubar)

        # --- Данные и состояния ---
        self.mm = None
        self.total = 0
        self.ecg = None
        self.r_peaks = np.array([])
        self.bad_regions = []        # список (start,end)
        self.region_patches = []     # Rectangle-объекты
        self.resizing = None         # (patch, 'left'/'right')
        self.notch_enabled = True
        self.manual_peak_mode = False  # режим ручного добавления/удаления пиков
        self.select_complex_mode = False  # режим выбора комплекса для отображения
        self.selected_idx = None
        self.manual_peaks = set()    # множество индексов ручно добавленных пиков
        
        # Переменные для кнопок
        self.manual_peak_button = None
        self.select_complex_button = None
        
        # Переменные для измерений
        self.measurement_mode = False
        self.measurement_points = []  # список точек измерений [(x, y, label), ...]
        self.measurement_lines = []   # список линий измерений
        
        # Переменные для автоматических маркеров ЭКГ
        self.ecg_markers = {}  # словарь с автоматическими маркерами {'P': (x, y), 'Q': (x, y), ...}
        self.manual_markers = {}  # словарь с ручно скорректированными маркерами
        self.marker_edit_mode = False  # режим редактирования маркеров
        
        # Переменные для сохранения данных усредненного комплекса
        self.t_avg_data = None
        self.avg_data = None
        
        # Кэш для ускорения отрисовки
        self._peak_cache = None
        self._cache_valid = False
        
        # Переменные для навигации по записи
        self.current_position = 0  # текущая позиция в секундах
        self.display_duration = 100  # длительность отображаемого фрагмента в секундах
        
        # Система кэширования и оптимизации
        self.cache = {}  # Кэш для обработанных фрагментов
        self.cache_lock = threading.Lock()
        self.processing_queue = queue.Queue()
        self.display_queue = queue.Queue()
        self.processing_thread = None
        self.last_update_time = 0
        self.update_interval = 0.1  # Обновление каждые 100мс
        self.executor = ThreadPoolExecutor(max_workers=2)
        self.is_processing = False
        
        # Флаг управления навигацией matplotlib
        self._navigation_disabled = False
        
        # Флаг для отслеживания первой отрисовки графика
        self._first_draw_complete = False
        
        # --- Параметры алгоритма Томсона ---
        self.sensitivity = tk.DoubleVar(value=1.0)
        self.min_rr_ms = tk.IntVar(value=200)
        self.window_ms = tk.IntVar(value=120)
        self.search_window_ms = tk.IntVar(value=50)
        self.auto_recompute = tk.BooleanVar(value=False)  # Контроль автоматического пересчета

        # --- Разметка grid ---
        self.rowconfigure(0, weight=3)  # Raw ECG
        self.rowconfigure(1, weight=0)  # Панель инструментов matplotlib
        self.rowconfigure(2, weight=2)  # Усредненный комплекс и спектр
        self.rowconfigure(3, weight=2)  # Выбранный комплекс
        self.rowconfigure(4, weight=1)  # Текстовое поле
        self.rowconfigure(5, weight=0)  # Панель кнопок
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        # --- Raw ECG с встроенной навигацией matplotlib ---
        self.fig_raw, self.ax_raw = plt.subplots(figsize=(12,3))
        plt.tight_layout()
        
        # Создаем фрейм для графика и панели инструментов
        graph_frame = ttk.Frame(self)
        graph_frame.grid(row=0, column=0, columnspan=2, sticky='nsew')
        graph_frame.rowconfigure(0, weight=1)
        graph_frame.columnconfigure(0, weight=1)
        
        # Включаем встроенные инструменты навигации matplotlib
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        self.canvas_raw = FigureCanvasTkAgg(self.fig_raw, master=graph_frame)
        self.canvas_raw.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        
        # Добавляем панель инструментов matplotlib в отдельный фрейм
        toolbar_frame = ttk.Frame(self)
        toolbar_frame.grid(row=1, column=0, columnspan=2, sticky='ew')
        toolbar_frame.columnconfigure(0, weight=1)
        
        self.toolbar_raw = NavigationToolbar2Tk(self.canvas_raw, toolbar_frame)
        self.toolbar_raw.grid(row=0, column=0, sticky='ew')
        
        # Метки поверх графика
        self.hr_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="ЧСС: -- bpm", background='#fff')
        self.hr_label.place(relx=0.01, rely=0.01)
        
        self.peaks_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="Пики: 0", background='#fff')
        self.peaks_label.place(relx=0.01, rely=0.05)
        
        self.progress_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="", background='#fff')
        self.progress_label.place(relx=0.01, rely=0.09)

        # Переменные для выделения региона
        self.region_selecting = False
        self.region_start = None
        self.region_end = None
        
        # Подключаем события мыши к canvas
        self.canvas_raw.mpl_connect('button_press_event', self.on_click)
        self.canvas_raw.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas_raw.mpl_connect('button_release_event', self.on_release)
        
        # Подключаем события навигации matplotlib для сброса режимов
        self.canvas_raw.mpl_connect('key_press_event', self.on_navigation_event)
        
        # Заменено на логику внутри основных обработчиков событий
        


        # --- Усреднённый комплекс ---
        self.fig_avg, self.ax_avg = plt.subplots(figsize=(4,2))
        plt.tight_layout()
        self.canvas_avg = FigureCanvasTkAgg(self.fig_avg, master=self)
        self.canvas_avg.get_tk_widget().grid(row=2, column=0, sticky='nsew', padx=5, pady=5)
        
        # Подключаем обработчик кликов для усредненного комплекса
        self.fig_avg.canvas.mpl_connect('button_press_event', self.on_avg_click)

        # --- Спектр ---
        self.fig_sp, self.ax_sp = plt.subplots(figsize=(4,2))
        plt.tight_layout()
        self.canvas_sp = FigureCanvasTkAgg(self.fig_sp, master=self)
        self.canvas_sp.get_tk_widget().grid(row=2, column=1, sticky='nsew', padx=5, pady=5)



        # --- Выбранный комплекс ---
        self.fig_sel, self.ax_sel = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_sel = FigureCanvasTkAgg(self.fig_sel, master=self)
        self.canvas_sel.get_tk_widget().grid(row=3, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

        # --- Текстовое поле для параметров ---
        text_frame = ttk.LabelFrame(self, text="Параметры и измерения")
        text_frame.grid(row=4, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        
        # Создаем текстовое поле с прокруткой
        text_container = ttk.Frame(text_frame)
        text_container.pack(fill='both', expand=True, padx=5, pady=5)
        
        self.text_widget = tk.Text(text_container, height=8, width=80, wrap='word')
        scrollbar = ttk.Scrollbar(text_container, orient='vertical', command=self.text_widget.yview)
        self.text_widget.configure(yscrollcommand=scrollbar.set)
        
        self.text_widget.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        # Кнопки для работы с текстом
        text_buttons = ttk.Frame(text_frame)
        text_buttons.pack(fill='x', padx=5, pady=2)
        
        ttk.Button(text_buttons, text="Очистить", command=self.clear_text).pack(side='left', padx=5)
        ttk.Button(text_buttons, text="Сохранить в файл", command=self.save_text_to_file).pack(side='left', padx=5)
        ttk.Button(text_buttons, text="Обновить параметры", command=self.update_text_parameters).pack(side='right', padx=5)

        
        
        # --- Панель кнопок ---
        ctrl = ttk.Frame(self)
        ctrl.grid(row=5, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        ttk.Button(ctrl, text='Notch 50 Hz', command=self.toggle_notch).pack(side='left', padx=5)
        self.region_selector_button = tk.Button(ctrl, text='Выделить регион', command=self.activate_region_selector,
                                               bg='lightgray', relief='raised', bd=2)
        self.region_selector_button.pack(side='left', padx=5)
        
        # Кнопки с визуальной индикацией состояния
        self.manual_peak_button = tk.Button(ctrl, text='Отметить комплекс', command=self.toggle_peak_selector,
                                           bg='lightgray', relief='raised', bd=2)
        self.manual_peak_button.pack(side='left', padx=5)
        
        self.select_complex_button = tk.Button(ctrl, text='Выбрать комплекс', command=self.toggle_select_complex,
                                               bg='lightgray', relief='raised', bd=2)
        self.select_complex_button.pack(side='left', padx=5)
        
        # Кнопка для измерений на усредненном комплексе
        self.measurement_button = tk.Button(ctrl, text='Измерения', command=self.toggle_measurement_mode,
                                           bg='lightgray', relief='raised', bd=2)
        self.measurement_button.pack(side='left', padx=5)
        
        # Кнопка для редактирования маркеров ЭКГ
        self.marker_edit_button = tk.Button(ctrl, text='Редактировать маркеры', command=self.toggle_marker_edit,
                                           bg='lightgray', relief='raised', bd=2)
        self.marker_edit_button.pack(side='left', padx=5)
        
        ttk.Button(ctrl, text='Обновить', command=self.compute).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Полная обработка', command=self.full_processing).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Copy Screenshot', command=self.copy_screenshot).pack(side='left', padx=5)

    # --- Методы GUI ---
    def open_file(self):
        fn = filedialog.askopenfilename(filetypes=[('BIN','*.bin')])
        if not fn:
            return
        try:
            self.mm, self.total = open_ecg_memmap(fn)
            self.clear_markup()
            
            # Загружаем полную запись для отображения в matplotlib
            self.load_full_recording()
            
        except Exception as e:
            messagebox.showerror('Ошибка', str(e))
            
    def load_full_recording(self):
        """Загружает полную запись для отображения в matplotlib"""
        if self.mm is None:
            return
            
        try:
            # Показываем прогресс
            self.progress_label.config(text="Загрузка записи...")
            self.update()
            
            # Извлекаем полную запись
            full_ecg = extract_channel(self.mm, self.total, 0, 0, self.total).astype(np.float64)
            full_ecg = full_ecg * ((2*2.4)/(2**24)) * 1e3
            
            # Применяем фильтры
            full_ecg = highpass_filter(full_ecg, FS, cutoff=0.5)
            if self.notch_enabled:
                full_ecg = notch_filter(full_ecg, FS)
            
            # Создаем временную шкалу
            t_full = np.arange(len(full_ecg)) / FS
            
            # Сохраняем данные для использования
            self.full_ecg = full_ecg
            self.t_full = t_full
            
            # Отрисовываем полную запись с пиками
            self.draw_raw()
            
            # Устанавливаем начальные границы просмотра - 100 секунд (только при первой загрузке)
            if not self._first_draw_complete:
                self.ax_raw.set_xlim(0, min(100, len(full_ecg) / FS))
                self.canvas_raw.draw()
            
            # Автоматически детектируем пики для отображения
            self.detect_peaks_for_display()
            
            # Очищаем прогресс
            self.progress_label.config(text="")
            
        except Exception as e:
            self.progress_label.config(text=f"Ошибка загрузки: {e}")
            print(f"Ошибка загрузки записи: {e}")
            
    def detect_peaks_for_display(self):
        """Детектирует пики для отображения на полной записи"""
        if not hasattr(self, 'full_ecg') or self.full_ecg is None:
            return
            
        try:
            self.progress_label.config(text="Детекция пиков...")
            self.update()
            
            # Детектируем пики на полной записи
            peaks, _, _, _, _ = detect_r_peaks(
                self.full_ecg, FS,
                sensitivity=self.sensitivity.get(),
                min_rr_ms=self.min_rr_ms.get(),
                window_ms=self.window_ms.get(),
                search_window_ms=self.search_window_ms.get(),
                notch_enabled=False  # Фильтры уже применены
            )
            
            # Сохраняем пики
            self.r_peaks = peaks
            
            # Обновляем отображение
            self.draw_raw()
            self.canvas_raw.draw()
            
            # Обновляем метки
            self.peaks_label.config(text=f"Пики: {len(self.r_peaks)}")
            
            # Вычисляем ЧСС только из годных пиков (исключаем плохие регионы)
            good_peaks = [r for r in self.r_peaks if not any(a<=r/FS<=b for a,b in self.bad_regions)]
            
            if len(good_peaks) > 1:
                rr_intervals = np.diff(np.array(good_peaks)) / FS * 1000  # в мс
                hr = 60000 / np.mean(rr_intervals)  # ударов в минуту
                self.hr_label.config(text=f"ЧСС: {hr:.1f} bpm")
            else:
                self.hr_label.config(text="ЧСС: -- bpm")
            
            self.progress_label.config(text="")
            
        except Exception as e:
            self.progress_label.config(text=f"Ошибка детекции: {e}")
            print(f"Ошибка детекции пиков: {e}")

    def toggle_notch(self):
        self.notch_enabled = not self.notch_enabled
        self.compute()

    def activate_region_selector(self):
        if self.mm is not None:
            # Переключаем режим выделения региона
            if not hasattr(self, 'region_selector_active'):
                self.region_selector_active = False
                
            self.region_selector_active = not self.region_selector_active
            
            # Отключаем другие режимы
            self.manual_peak_mode = False
            self.select_complex_mode = False
            self.measurement_mode = False
            self.marker_edit_mode = False
            
            # Обновляем состояние кнопок
            if self.region_selector_active:
                self.region_selector_button.config(bg='lightcoral', relief='sunken')
                self.manual_peak_button.config(bg='lightgray', relief='raised')
                self.select_complex_button.config(bg='lightgray', relief='raised')
                self.measurement_button.config(bg='lightgray', relief='raised')
                self.marker_edit_button.config(bg='lightgray', relief='raised')
                
                # Полностью отключаем инструменты навигации matplotlib
                self.disable_matplotlib_navigation()
            else:
                self.region_selector_button.config(bg='lightgray', relief='raised')
                
                # Включаем инструменты навигации matplotlib
                self.enable_matplotlib_navigation()

    def toggle_peak_selector(self):
        """Переключение режима ручного добавления/удаления пиков"""
        self.manual_peak_mode = not self.manual_peak_mode
        self.select_complex_mode = False  # отключаем другой режим
        self.measurement_mode = False
        self.marker_edit_mode = False
        
        # Обновляем состояние кнопок
        if self.manual_peak_mode:
            self.manual_peak_button.config(bg='lightgreen', relief='sunken')
            self.region_selector_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            self.measurement_button.config(bg='lightgray', relief='raised')
            self.marker_edit_button.config(bg='lightgray', relief='raised')
        else:
            self.manual_peak_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector при активном режиме ручной разметки
        self.span.set_active(False)
        
        # Полностью отключаем инструменты навигации matplotlib
        self.disable_matplotlib_navigation()

    def toggle_select_complex(self):
        """Переключение режима выбора комплекса для отображения"""
        self.select_complex_mode = not self.select_complex_mode
        self.manual_peak_mode = False  # отключаем другой режим
        self.measurement_mode = False
        self.marker_edit_mode = False
        
        # Обновляем состояние кнопок
        if self.select_complex_mode:
            self.select_complex_button.config(bg='lightblue', relief='sunken')
            self.region_selector_button.config(bg='lightgray', relief='raised')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.measurement_button.config(bg='lightgray', relief='raised')
            self.marker_edit_button.config(bg='lightgray', relief='raised')
        else:
            self.select_complex_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector при активном режиме выбора комплекса
        self.span.set_active(False)
        
        # Полностью отключаем инструменты навигации matplotlib
        self.disable_matplotlib_navigation()

    def toggle_measurement_mode(self):
        """Переключение режима измерений на усредненном комплексе"""
        self.measurement_mode = not self.measurement_mode
        self.manual_peak_mode = False
        self.select_complex_mode = False
        self.marker_edit_mode = False
        
        # Обновляем состояние кнопок
        if self.measurement_mode:
            self.measurement_button.config(bg='lightyellow', relief='sunken')
            self.region_selector_button.config(bg='lightgray', relief='raised')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            self.marker_edit_button.config(bg='lightgray', relief='raised')
            messagebox.showinfo("Режим измерений", 
                              "Включен режим измерений.\n"
                              "Левый клик - добавить маркер ЭКГ (P_начало, P_конец, Q, R, S, T_начало, T_конец)\n"
                              "Правый клик - удалить ближайший маркер")
        else:
            self.measurement_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector
        self.span.set_active(False)
        
        # Отключаем инструменты навигации matplotlib
        self.toolbar_raw.pan()
        self.toolbar_raw.zoom()

    def toggle_marker_edit(self):
        """Переключение режима редактирования маркеров ЭКГ"""
        self.marker_edit_mode = not self.marker_edit_mode
        self.measurement_mode = False
        self.manual_peak_mode = False
        self.select_complex_mode = False
        
        # Обновляем состояние кнопок
        if self.marker_edit_mode:
            self.marker_edit_button.config(bg='lightcoral', relief='sunken')
            self.region_selector_button.config(bg='lightgray', relief='raised')
            self.measurement_button.config(bg='lightgray', relief='raised')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            messagebox.showinfo("Редактирование маркеров", 
                              "Включен режим редактирования маркеров ЭКГ.\n"
                              "Левый клик - переместить маркер\n"
                              "Правый клик - удалить маркер")
        else:
            self.marker_edit_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector
        self.span.set_active(False)
        
        # Полностью отключаем инструменты навигации matplotlib
        self.disable_matplotlib_navigation()
        
    def on_navigation_event(self, event):
        """Обработчик событий навигации matplotlib - сбрасывает активные режимы"""
        # Если пользователь использует инструменты навигации matplotlib, сбрасываем наши режимы
        if hasattr(event, 'key') and event.key in ['p', 'o', 'home', 'backspace']:
            self.reset_all_modes()
            

        

                
    def reset_all_modes(self):
        """Сбрасывает все активные режимы"""
        self.manual_peak_mode = False
        self.select_complex_mode = False
        self.measurement_mode = False
        self.marker_edit_mode = False
        self.region_selector_active = False
        
        # Обновляем состояние кнопок
        self.region_selector_button.config(bg='lightgray', relief='raised')
        self.manual_peak_button.config(bg='lightgray', relief='raised')
        self.select_complex_button.config(bg='lightgray', relief='raised')
        self.measurement_button.config(bg='lightgray', relief='raised')
        self.marker_edit_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector
        self.span.set_active(False)
        
        # Включаем инструменты навигации matplotlib
        self.enable_matplotlib_navigation()
        
    def disable_matplotlib_navigation(self):
        """Полностью отключает инструменты навигации matplotlib"""
        # Устанавливаем флаг, что навигация отключена
        self._navigation_disabled = True
        
        # Принудительно сбрасываем все активные режимы навигации
        try:
            # Отключаем режим панорамирования если активен
            if hasattr(self.toolbar_raw, '_active') and self.toolbar_raw._active == 'PAN':
                self.toolbar_raw.pan()
            
            # Отключаем режим зума если активен  
            if hasattr(self.toolbar_raw, '_active') and self.toolbar_raw._active == 'ZOOM':
                self.toolbar_raw.zoom()
        except:
            pass
        
        # Блокируем кнопки навигации в toolbar
        for child in self.toolbar_raw.winfo_children():
            if hasattr(child, 'config'):
                try:
                    child.config(state='disabled')
                except:
                    pass
        
    def enable_matplotlib_navigation(self):
        """Включает инструменты навигации matplotlib"""
        # Сбрасываем флаг отключения навигации
        self._navigation_disabled = False
        
        # Разблокируем кнопки навигации в toolbar
        for child in self.toolbar_raw.winfo_children():
            if hasattr(child, 'config'):
                try:
                    child.config(state='normal')
                except:
                    pass

    def clear_text(self):
        """Очистить текстовое поле"""
        self.text_widget.delete(1.0, tk.END)

    def save_text_to_file(self):
        """Сохранить содержимое текстового поля в файл"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Сохранить параметры"
        )
        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(self.text_widget.get(1.0, tk.END))
                messagebox.showinfo("Успех", "Файл сохранен успешно!")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить файл: {str(e)}")

    def update_text_parameters(self):
        """Обновить параметры в текстовом поле"""
        if not hasattr(self, 'full_ecg') or self.full_ecg is None or len(self.r_peaks) == 0:
            messagebox.showwarning("Предупреждение", "Нет данных для анализа")
            return
        
        # Вычисляем параметры
        good_peaks = [r for r in self.r_peaks if not any(a<=r/FS<=b for a,b in self.bad_regions)]
        
        # Вычисляем общую длительность записи
        total_duration = self.total / FS  # в секундах
        
        if len(good_peaks) > 1:
            rr_intervals = np.diff(np.array(good_peaks)/FS)  # в секундах
            hr_values = 60 / rr_intervals  # ЧСС в ударах в минуту
            
            hr_mean = np.mean(hr_values)
            hr_min = np.min(hr_values)
            hr_max = np.max(hr_values)
            hr_std = np.std(hr_values)
            
            # Формируем текст
            text_content = f"""=== АНАЛИЗ ЭКГ ===
Дата: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

ОБЩИЕ ПАРАМЕТРЫ:
• Длительность записи: {total_duration:.1f} сек ({total_duration/60:.1f} мин)
• Количество обнаруженных пиков: {len(self.r_peaks)}
• Количество годных пиков: {len(good_peaks)}
• Количество ручных пиков: {len(self.manual_peaks)}

ЧАСТОТА СЕРДЕЧНЫХ СОКРАЩЕНИЙ (ЧСС):
• Средняя ЧСС: {hr_mean:.1f} ± {hr_std:.1f} уд/мин
• Минимальная ЧСС: {hr_min:.1f} уд/мин
• Максимальная ЧСС: {hr_max:.1f} уд/мин
• Вариабельность ЧСС: {hr_std:.1f} уд/мин

RR ИНТЕРВАЛЫ:
• Средний RR интервал: {np.mean(rr_intervals)*1000:.1f} мс
• Минимальный RR интервал: {np.min(rr_intervals)*1000:.1f} мс
• Максимальный RR интервал: {np.max(rr_intervals)*1000:.1f} мс
• Стандартное отклонение RR: {np.std(rr_intervals)*1000:.1f} мс

ПАРАМЕТРЫ АЛГОРИТМА ТОМСОНА:
• Чувствительность: {self.sensitivity.get():.2f}
• Минимальное расстояние между пиками: {self.min_rr_ms.get()} мс
• Окно интегрирования: {self.window_ms.get()} мс
• Окно поиска: {self.search_window_ms.get()} мс

ИЗМЕРЕНИЯ НА УСРЕДНЕННОМ КОМПЛЕКСЕ:
"""
            
            # Добавляем маркеры ЭКГ
            all_markers = {**self.ecg_markers, **self.manual_markers}
            if all_markers:
                text_content += "\nМАРКЕРЫ ЭКГ:\n"
                marker_order = ['P_начало', 'P_конец', 'Q', 'R', 'S', 'T_начало', 'T_конец']
                
                # Определяем точку отсчета времени (P_начало = 0 мс)
                p_start_time = None
                if 'P_начало' in all_markers:
                    p_start_time = all_markers['P_начало'][0]  # время P_начало в мс относительно центра (R-пик)
                
                for marker_name in marker_order:
                    if marker_name in all_markers:
                        x, y = all_markers[marker_name]  # x в мс относительно центра комплекса (R-пик)
                        status = "(ручная)" if marker_name in self.manual_markers else "(авто)"
                        
                        # Вычисляем абсолютное время от P_начало
                        if p_start_time is not None:
                            abs_time = x - p_start_time  # уже в мс от P_начало
                            text_content += f"• {marker_name} {status}: {abs_time:.1f} мс, {y:.2f} мВ\n"
                        else:
                            # Если P_начало не установлен, показываем относительно центра комплекса (R-пик)
                            text_content += f"• {marker_name} {status}: {x:.1f} мс (от R-пика), {y:.2f} мВ\n"
                
                # Вычисляем интервалы и абсолютные времена
                text_content += "\nИНТЕРВАЛЫ И АБСОЛЮТНЫЕ ВРЕМЕНА (от P_начало = 0 мс):\n"
                
                # Длительность P волны
                if 'P_начало' in all_markers and 'P_конец' in all_markers:
                    p_duration = all_markers['P_конец'][0] - all_markers['P_начало'][0]  # уже в мс
                    text_content += f"• Длительность P волны: {p_duration:.1f} мс\n"
                
                # PQ интервал (от конца P до начала Q)
                if 'P_конец' in all_markers and 'Q' in all_markers:
                    pq_interval = all_markers['Q'][0] - all_markers['P_конец'][0]  # уже в мс
                    text_content += f"• PQ интервал (изоэлектрический): {pq_interval:.1f} мс\n"
                
                # QRS длительность
                if 'Q' in all_markers and 'S' in all_markers:
                    qrs_duration = all_markers['S'][0] - all_markers['Q'][0]  # уже в мс
                    text_content += f"• QRS длительность: {qrs_duration:.1f} мс\n"
                
                # ST интервал (от S до начала T)
                if 'S' in all_markers and 'T_начало' in all_markers:
                    st_interval = all_markers['T_начало'][0] - all_markers['S'][0]  # уже в мс
                    text_content += f"• ST интервал (изоэлектрический): {st_interval:.1f} мс\n"
                
                # Длительность T волны
                if 'T_начало' in all_markers and 'T_конец' in all_markers:
                    t_duration = all_markers['T_конец'][0] - all_markers['T_начало'][0]  # уже в мс
                    text_content += f"• Длительность T волны: {t_duration:.1f} мс\n"
                
                # PR интервал (от начала P до R)
                if 'P_начало' in all_markers and 'R' in all_markers:
                    pr_interval = all_markers['R'][0] - all_markers['P_начало'][0]  # уже в мс
                    text_content += f"• PR интервал (общий): {pr_interval:.1f} мс\n"
                
                # QT интервал (от Q до конца T)
                if 'Q' in all_markers and 'T_конец' in all_markers:
                    qt_interval = all_markers['T_конец'][0] - all_markers['Q'][0]  # уже в мс
                    text_content += f"• QT интервал (общий): {qt_interval:.1f} мс\n"
                
                # Полная длительность кардиоцикла (от P_начало до T_конец)
                if 'P_начало' in all_markers and 'T_конец' in all_markers:
                    total_duration = all_markers['T_конец'][0] - all_markers['P_начало'][0]  # уже в мс
                    text_content += f"• Полная длительность кардиоцикла: {total_duration:.1f} мс\n"
            
            # Добавляем измерения
            if self.measurement_points:
                text_content += "\nДОПОЛНИТЕЛЬНЫЕ ИЗМЕРЕНИЯ:\n"
                for i, (x, y, label) in enumerate(self.measurement_points, 1):
                    text_content += f"• {label}: {x:.1f} мс, {y:.2f} мВ\n"
            
            # Добавляем в текстовое поле
            self.text_widget.delete(1.0, tk.END)
            self.text_widget.insert(1.0, text_content)
        else:
            messagebox.showwarning("Предупреждение", "Недостаточно пиков для анализа")

    def on_avg_click(self, event):
        """Обработчик кликов на усредненном комплексе"""
        if event.inaxes != self.ax_avg:
            return
        
        # Режим измерений
        if self.measurement_mode:
            if event.button == 1:  # Левый клик - добавить маркер
                # Определяем какой маркер добавить (P_начало, P_конец, Q, R, S, T_начало, T_конец)
                existing_markers = set(self.manual_markers.keys())
                available_markers = ['P_начало', 'P_конец', 'Q', 'R', 'S', 'T_начало', 'T_конец']
                next_marker = None
                
                for marker in available_markers:
                    if marker not in existing_markers:
                        next_marker = marker
                        break
                
                if next_marker:
                    # Добавляем маркер ЭКГ
                    self.manual_markers[next_marker] = (event.xdata, event.ydata)
                    
                    # Перерисовываем график
                    self.draw_average_complex()
                
            elif event.button == 3:  # Правый клик - удалить ближайший маркер
                all_markers = {**self.ecg_markers, **self.manual_markers}
                if all_markers:
                    distances = [(abs(event.xdata - x), marker_name) for marker_name, (x, y) in all_markers.items()]
                    distances.sort()
                    nearest_marker = distances[0][1]
                    
                    # Удаляем маркер
                    if nearest_marker in self.manual_markers:
                        del self.manual_markers[nearest_marker]
                    if nearest_marker in self.ecg_markers:
                        del self.ecg_markers[nearest_marker]
                    
                    # Перерисовываем график
                    self.draw_average_complex()
        
        # Режим редактирования маркеров ЭКГ
        elif self.marker_edit_mode:
            if event.button == 1:  # Левый клик - переместить маркер
                # Находим ближайший маркер ЭКГ
                all_markers = {**self.ecg_markers, **self.manual_markers}
                if all_markers:
                    distances = [(abs(event.xdata - x), marker_name) for marker_name, (x, y) in all_markers.items()]
                    distances.sort()
                    nearest_marker = distances[0][1]
                    
                    # Перемещаем ближайший маркер
                    self.manual_markers[nearest_marker] = (event.xdata, event.ydata)
                    
                    # Перерисовываем график
                    self.draw_average_complex()
                    
            elif event.button == 3:  # Правый клик - удалить маркер
                all_markers = {**self.ecg_markers, **self.manual_markers}
                if all_markers:
                    distances = [(abs(event.xdata - x), marker_name) for marker_name, (x, y) in all_markers.items()]
                    distances.sort()
                    nearest_marker = distances[0][1]
                    
                    # Удаляем маркер
                    if nearest_marker in self.manual_markers:
                        del self.manual_markers[nearest_marker]
                    if nearest_marker in self.ecg_markers:
                        del self.ecg_markers[nearest_marker]
                    
                    # Перерисовываем график
                    self.draw_average_complex()

    def draw_average_complex(self, t_avg=None, avg=None):
        """Отрисовка усредненного комплекса с маркерами"""
        self.ax_avg.cla()
        
        # Сохраняем данные если они переданы
        if t_avg is not None and avg is not None:
            self.t_avg_data = t_avg
            self.avg_data = avg
        
        # Используем сохраненные данные если новые не переданы
        if self.t_avg_data is not None and self.avg_data is not None:
            self.ax_avg.plot(self.t_avg_data, self.avg_data, linewidth=1, color='blue')
            
            # Отображаем маркеры ЭКГ
            self.draw_ecg_markers(self.t_avg_data, self.avg_data)
        
        self.ax_avg.set_title('Усреднённый комплекс')
        self.ax_avg.set_xlabel('Время (мс)')
        self.ax_avg.set_ylabel('Амплитуда (мВ)')
        self.ax_avg.grid(True, alpha=0.3)
        self.canvas_avg.draw()

    def draw_ecg_markers(self, t_avg, avg):
        """Отрисовка маркеров ЭКГ"""
        # Цвета для разных маркеров
        colors = {
            'P_начало': 'green', 'P_конец': 'darkgreen', 
            'Q': 'orange', 'R': 'red', 'S': 'purple', 
            'T_начало': 'blue', 'T_конец': 'darkblue'
        }
        
        # Объединяем автоматические и ручные маркеры
        all_markers = {**self.ecg_markers, **self.manual_markers}
        
        for marker_name, (x, y) in all_markers.items():
            if marker_name in colors:
                color = colors[marker_name]
                # Ручные маркеры рисуем больше
                size = 8 if marker_name in self.manual_markers else 6
                self.ax_avg.plot(x, y, 'o', color=color, markersize=size, markeredgecolor='black', markeredgewidth=1)
                self.ax_avg.annotate(marker_name, (x, y), 
                                   xytext=(5, 5), textcoords='offset points',
                                   fontsize=10, fontweight='bold', color=color,
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    def on_param_change(self, event=None):
        """Обработчик изменения параметров алгоритма Томсона"""
        # Обновляем метки
        self.sensitivity_label.config(text=f"{self.sensitivity.get():.1f}")
        self.min_rr_label.config(text=f"{self.min_rr_ms.get()}")
        self.window_label.config(text=f"{self.window_ms.get()}")
        self.search_label.config(text=f"{self.search_window_ms.get()}")
        
        # Пересчитываем только если включен автоматический режим
        if self.auto_recompute.get() and self.mm is not None:
            self.compute()

    def reset_tompson_params(self):
        """Сброс параметров алгоритма Томсона к значениям по умолчанию"""
        self.sensitivity.set(1.0)
        self.min_rr_ms.set(200)
        self.window_ms.set(120)
        self.search_window_ms.set(50)
        
        # Обновляем метки
        self.sensitivity_label.config(text="1.0")
        self.min_rr_label.config(text="200")
        self.window_label.config(text="120")
        self.search_label.config(text="50")
        
        # Пересчитываем
        if self.mm is not None:
            self.compute()

    def save_params(self):
        """Сохранение параметров алгоритма Томсона в файл"""
        import json
        
        params = {
            'sensitivity': self.sensitivity.get(),
            'min_rr_ms': self.min_rr_ms.get(),
            'window_ms': self.window_ms.get(),
            'search_window_ms': self.search_window_ms.get()
        }
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Сохранить параметры алгоритма Томсона"
        )
        
        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    json.dump(params, f, indent=2, ensure_ascii=False)
                messagebox.showinfo("Успех", "Параметры сохранены успешно!")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить параметры: {str(e)}")

    def load_params(self):
        """Загрузка параметров алгоритма Томсона из файла"""
        import json
        
        filename = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Загрузить параметры алгоритма Томсона"
        )
        
        if filename:
            try:
                with open(filename, 'r', encoding='utf-8') as f:
                    params = json.load(f)
                
                # Устанавливаем параметры
                self.sensitivity.set(params.get('sensitivity', 1.0))
                self.min_rr_ms.set(params.get('min_rr_ms', 200))
                self.window_ms.set(params.get('window_ms', 120))
                self.search_window_ms.set(params.get('search_window_ms', 50))
                
                # Обновляем метки
                self.sensitivity_label.config(text=f"{self.sensitivity.get():.1f}")
                self.min_rr_label.config(text=f"{self.min_rr_ms.get()}")
                self.window_label.config(text=f"{self.window_ms.get()}")
                self.search_label.config(text=f"{self.search_window_ms.get()}")
                
                # Пересчитываем
                if self.mm is not None:
                    self.compute()
                    
                messagebox.showinfo("Успех", "Параметры загружены успешно!")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось загрузить параметры: {str(e)}")

    def full_processing(self):
        """Полная обработка файла с последовательной детекцией пиков"""
        if self.mm is None:
            messagebox.showwarning("Предупреждение", "Сначала загрузите файл ЭКГ")
            return
        
        try:
            # Показываем диалог прогресса
            progress_window = tk.Toplevel(self)
            progress_window.title("Обработка файла")
            progress_window.geometry("300x150")
            progress_window.transient(self)
            progress_window.grab_set()
            
            progress_label = ttk.Label(progress_window, text="Обработка файла...")
            progress_label.pack(pady=20)
            
            progress_bar = ttk.Progressbar(progress_window, mode='determinate', maximum=100)
            progress_bar.pack(pady=10, padx=20, fill='x')
            
            def update_progress(progress):
                progress_bar['value'] = progress
                progress_label.config(text=f"Обработано: {progress:.1f}%")
                progress_window.update()
            
            # Обновляем GUI
            progress_window.update()
            
            # Выполняем последовательную обработку
            self.r_peaks = detect_r_peaks_streaming(
                self.mm, self.total, 0, FS,
                chunk_size=10000,  # 10 секунд
                overlap=2000,      # 2 секунды перекрытия
                progress_callback=update_progress,
                sensitivity=self.sensitivity.get(),
                min_rr_ms=self.min_rr_ms.get(),
                window_ms=self.window_ms.get(),
                search_window_ms=self.search_window_ms.get(),
                notch_enabled=self.notch_enabled
            )
            
            # Закрываем окно прогресса
            progress_window.destroy()
            
            # Обновляем отображение
            self.compute()
            
            messagebox.showinfo("Успех", f"Обработано {len(self.r_peaks)} пиков")
            
        except Exception as e:
                            messagebox.showerror("Ошибка", f"Ошибка при обработке: {str(e)}")

    def show_tompson_settings(self):
        """Показать окно настроек алгоритма Томсона"""
        settings_window = tk.Toplevel(self)
        settings_window.title("Параметры алгоритма Томсона")
        settings_window.geometry("600x300")
        settings_window.transient(self)
        settings_window.grab_set()
        
        # Основной фрейм
        main_frame = ttk.Frame(settings_window, padding="10")
        main_frame.pack(fill='both', expand=True)
        
        # Заголовок
        ttk.Label(main_frame, text="Настройки алгоритма Томсона", font=('Arial', 12, 'bold')).pack(pady=(0, 20))
        
        # Первая строка параметров
        param_row1 = ttk.Frame(main_frame)
        param_row1.pack(fill='x', pady=5)
        
        ttk.Label(param_row1, text="Чувствительность:", width=20).pack(side='left', padx=5)
        sensitivity_scale = ttk.Scale(param_row1, from_=0.1, to=2.0, variable=self.sensitivity, 
                                     orient='horizontal', length=200, command=self.on_param_change)
        sensitivity_scale.pack(side='left', padx=5)
        self.sensitivity_label = ttk.Label(param_row1, text="1.0", width=8)
        self.sensitivity_label.pack(side='left', padx=5)
        
        ttk.Label(param_row1, text="Мин. RR (мс):", width=15).pack(side='left', padx=5)
        min_rr_scale = ttk.Scale(param_row1, from_=100, to=500, variable=self.min_rr_ms, 
                                orient='horizontal', length=150, command=self.on_param_change)
        min_rr_scale.pack(side='left', padx=5)
        self.min_rr_label = ttk.Label(param_row1, text="200", width=8)
        self.min_rr_label.pack(side='left', padx=5)
        
        # Вторая строка параметров
        param_row2 = ttk.Frame(main_frame)
        param_row2.pack(fill='x', pady=5)
        
        ttk.Label(param_row2, text="Окно интегрирования (мс):", width=20).pack(side='left', padx=5)
        window_scale = ttk.Scale(param_row2, from_=50, to=200, variable=self.window_ms, 
                                orient='horizontal', length=200, command=self.on_param_change)
        window_scale.pack(side='left', padx=5)
        self.window_label = ttk.Label(param_row2, text="120", width=8)
        self.window_label.pack(side='left', padx=5)
        
        ttk.Label(param_row2, text="Окно поиска (мс):", width=15).pack(side='left', padx=5)
        search_scale = ttk.Scale(param_row2, from_=20, to=100, variable=self.search_window_ms, 
                                orient='horizontal', length=150, command=self.on_param_change)
        search_scale.pack(side='left', padx=5)
        self.search_label = ttk.Label(param_row2, text="50", width=8)
        self.search_label.pack(side='left', padx=5)
        
        # Чекбокс для автоматического пересчета
        auto_frame = ttk.Frame(main_frame)
        auto_frame.pack(fill='x', pady=10)
        ttk.Checkbutton(auto_frame, text="Автоматический пересчет при изменении параметров", 
                       variable=self.auto_recompute).pack(side='left', padx=5)
        
        # Кнопки
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill='x', pady=20)
        
        ttk.Button(button_frame, text="Применить", command=lambda: [self.compute(), settings_window.destroy()]).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Отмена", command=settings_window.destroy).pack(side='left', padx=5)
        ttk.Button(button_frame, text="Сброс к значениям по умолчанию", command=self.reset_tompson_params).pack(side='right', padx=5)
        
        # Описание параметров
        desc_frame = ttk.LabelFrame(main_frame, text="Описание параметров", padding="10")
        desc_frame.pack(fill='x', pady=10)
        
        desc_text = """
• Чувствительность (0.1-2.0): Основной параметр контроля детекции. Меньшие значения = более строгая детекция.
• Мин. RR (100-500 мс): Минимальный интервал между пиками. Помогает избежать множественных детекций.
• Окно интегрирования (50-200 мс): Размер окна для сглаживания дифференцированного сигнала.
• Окно поиска (20-100 мс): Окно для поиска точного положения R-пика в исходном сигнале.
        """
        ttk.Label(desc_frame, text=desc_text, justify='left').pack(anchor='w')

    def clear_markup(self):
        # удалить патчи
        for p in self.region_patches:
            try: p.remove()
            except: pass
        self.region_patches = []
        self.bad_regions = []
        self.manual_peaks = set()
        self.selected_idx = None
        self.measurement_points = []
        self.ecg_markers = {}
        self.manual_markers = {}
        self.t_avg_data = None
        self.avg_data = None
        
        # Сбрасываем флаг первой отрисовки для нового файла
        self._first_draw_complete = False
        
        self.compute()



    def on_click(self, ev):
        if ev.inaxes != self.ax_raw:
            return
            
        # Обработка выделения региона
        if getattr(self, 'region_selector_active', False):
            if ev.button == 1:  # Левый клик
                if not self.region_selecting:
                    # Начинаем выделение региона
                    self.region_selecting = True
                    self.region_start = ev.xdata
                    print(f"Начало выделения региона: {self.region_start:.2f} сек")
                else:
                    # Завершаем выделение региона
                    self.region_selecting = False
                    self.region_end = ev.xdata
                    print(f"Конец выделения региона: {self.region_end:.2f} сек")
                    
                    # Создаем регион
                    start, end = min(self.region_start, self.region_end), max(self.region_start, self.region_end)
                    self.bad_regions.append((start, end))
                    
                    # Создаем визуальный патч
                    if hasattr(self, 'full_ecg') and self.full_ecg is not None:
                        y_min = np.min(self.full_ecg)
                        y_max = np.max(self.full_ecg)
                    else:
                        y_min, y_max = -1, 1
                    
                    patch = Rectangle((start, y_min), end-start, y_max-y_min,
                                      color='red', alpha=0.3)
                    self.ax_raw.add_patch(patch)
                    self.region_patches.append(patch)
                    self.canvas_raw.draw()
                    
                    # Сбрасываем режим выделения
                    self.region_selector_active = False
                    self.region_selector_button.config(bg='lightgray', relief='raised')
                    self.enable_matplotlib_navigation()
            
            # ВАЖНО: Блокируем дальнейшую обработку события matplotlib при активном режиме выделения региона
            return
            
        # попытка редактирования region
        for idx, patch in enumerate(self.region_patches):
            contains, _ = patch.contains(ev)
            if contains:
                x0, w = patch.get_x(), patch.get_width()
                tol = w * 0.05
                if abs(ev.xdata - x0) < tol:
                    self.resizing = (patch, 'left', idx)
                    return
                elif abs(ev.xdata - (x0 + w)) < tol:
                    self.resizing = (patch, 'right', idx)
                    return
                else:
                    patch.remove()
                    del self.region_patches[idx]
                    del self.bad_regions[idx]
                    self.canvas_raw.draw()
                    return
        
        # Ручное добавление/удаление пиков
        if self.manual_peak_mode:
            t = ev.xdata  # Это уже абсолютное время
            # Конвертируем абсолютное время в индекс сэмпла
            sample_idx = int(t * FS)
            
            # Проверяем границы
            if sample_idx < 0 or sample_idx >= len(self.full_ecg):
                return
            
            if ev.button == 3:  # Правый клик - добавить пик
                # Ищем локальный максимум в окрестности клика
                search_window = int(0.1 * FS)  # 100 мс окно поиска
                start = max(0, sample_idx - search_window)
                end = min(len(self.full_ecg), sample_idx + search_window)
                
                # Находим локальный максимум
                local_max_idx = start + np.argmax(self.full_ecg[start:end])
                
                # Добавляем пик в список
                self.r_peaks = np.append(self.r_peaks, local_max_idx)
                self.r_peaks = np.sort(self.r_peaks)  # сортируем
                # Находим новый индекс добавленного пика
                new_idx = np.where(self.r_peaks == local_max_idx)[0][0]
                self.manual_peaks.add(new_idx)  # запоминаем как ручной
                
                self._invalidate_cache()
                self.draw_raw()
                
            elif ev.button == 1:  # Левый клик - удалить пик
                # Ищем ближайший пик
                if len(self.r_peaks) > 0:
                    # Конвертируем время клика в полную частоту дискретизации
                    click_time_full = sample_idx / FS
                    distances = np.abs(self.r_peaks/FS - click_time_full)
                    nearest_idx = np.argmin(distances)
                    
                    if distances[nearest_idx] < 0.1:  # 100 мс допуск
                        # Удаляем пик
                        self.r_peaks = np.delete(self.r_peaks, nearest_idx)
                        # Обновляем индексы ручных пиков
                        self.manual_peaks = {i-1 if i > nearest_idx else i for i in self.manual_peaks}
                        self.manual_peaks.discard(nearest_idx)
                        
                        # Обновляем индекс выбранного пика
                        if self.selected_idx is not None:
                            if self.selected_idx == nearest_idx:
                                self.selected_idx = None
                            elif self.selected_idx > nearest_idx:
                                self.selected_idx -= 1
                        
                        self._invalidate_cache()
                        self.draw_raw()
        
        # Выбор комплекса для отображения
        elif self.select_complex_mode and ev.button == 1:
            if len(self.r_peaks) > 0:
                # ev.xdata уже содержит абсолютное время
                click_time_full = ev.xdata
                
                # Ищем ближайший пик
                distances = np.abs(self.r_peaks/FS - click_time_full)
                nearest_idx = np.argmin(distances)
                
                if distances[nearest_idx] < 0.1:  # 100 мс допуск
                    self.selected_idx = nearest_idx
                    self.draw_selected()
                    self.draw_raw()  # Обновляем отображение для показа выбранного пика
                    
        # Режим измерений на усредненном комплексе
        elif self.measurement_mode:
            if ev.button == 1:  # Левый клик - добавить маркер
                # Здесь должна быть логика для добавления маркеров на усредненном комплексе
                pass
            elif ev.button == 3:  # Правый клик - удалить маркер
                # Здесь должна быть логика для удаления маркеров
                pass
                
        # Режим редактирования маркеров
        elif self.marker_edit_mode:
            if ev.button == 1:  # Левый клик - добавить/переместить маркер
                # Здесь должна быть логика для редактирования маркеров
                pass
            elif ev.button == 3:  # Правый клик - удалить маркер
                # Здесь должна быть логика для удаления маркеров
                pass

    def on_motion(self, ev):
        # Блокируем события matplotlib при активном режиме выделения региона
        if getattr(self, 'region_selector_active', False):
            return
            
        if not self.resizing or ev.inaxes!=self.ax_raw:
            return
        patch, side, idx = self.resizing
        x0, w = patch.get_x(), patch.get_width()
        if side=='left':
            new_x0 = min(ev.xdata, x0+w-0.01)
            patch.set_x(new_x0)
            patch.set_width(x0+w - new_x0)
            self.bad_regions[idx] = (new_x0, new_x0 + patch.get_width())
        else:
            new_w = max(ev.xdata - x0, 0.01)
            patch.set_width(new_w)
            self.bad_regions[idx] = (x0, x0 + new_w)
        self.canvas_raw.draw()

    def on_release(self, ev):
        # Блокируем события matplotlib при активном режиме выделения региона
        if getattr(self, 'region_selector_active', False):
            return
            
        self.resizing = None

    def compute(self):
        if self.mm is None:
            return
        
        # Показываем прогресс
        self.progress_label.config(text="Обработка всей записи...")
        self.update()
        
        # Обрабатываем всю запись последовательно
        self.r_peaks = detect_r_peaks_streaming(
            self.mm, self.total, 0, FS,
            chunk_size=CHUNK_SIZE_MS,
            overlap=OVERLAP_MS,
            progress_callback=self.update_progress,
            sensitivity=self.sensitivity.get(),
            min_rr_ms=self.min_rr_ms.get(),
            window_ms=self.window_ms.get(),
            search_window_ms=self.search_window_ms.get(),
            notch_enabled=self.notch_enabled
        )
        
        # HR по годным пикам (исключаем плохие регионы)
        good = [r for r in self.r_peaks if not any(a<=r/FS<=b for a,b in self.bad_regions)]
        
        # Обновляем информацию о пиках
        self.peaks_label.config(text=f"Пики: {len(self.r_peaks)}")
        
        if len(good)>1:
            rr = np.diff(np.array(good)/FS)
            hr = (60/rr).mean()
            self.hr_label.config(text=f"ЧСС: {hr:.0f} bpm")
        else:
            self.hr_label.config(text="ЧСС: -- bpm")

        # Обновляем отображение пиков на графике
        self.draw_raw()
        
        # Создаем усредненный комплекс из всех пиков
        t_avg, avg = average_complex_from_peaks(self.mm, self.total, good, FS)
        
        # Создаем маркеры ЭКГ для ручного позиционирования
        if avg is not None:
            self.ecg_markers = create_ecg_markers(t_avg, avg)
        
        self.draw_average_complex(t_avg, avg)

        self.ax_sp.cla()
        if avg is not None:
            ac = avg - np.mean(avg)
            w = np.hanning(len(ac))
            x = ac * w
            yf = np.fft.rfft(x)
            xf = np.fft.rfftfreq(len(x), 1/FS)
            amp = np.abs(yf) / len(x)
            mask = xf <= 150
            self.ax_sp.plot(xf[mask], amp[mask], linewidth=1)
        self.ax_sp.set_title('Спектр усреднённого комплекса')
        self.canvas_sp.draw()

        self.draw_selected()
        
        # Очищаем прогресс
        self.progress_label.config(text="")
        
        # Инвалидируем кэш
        self._cache_valid = False
    
    def update_progress(self, progress):
        """Обновление индикатора прогресса"""
        self.progress_label.config(text=f"Прогресс: {progress:.1f}%")
        self.update()
    
    def _invalidate_cache(self):
        """Инвалидация кэша при изменении данных"""
        self._cache_valid = False
        
    def get_cache_key(self, position, duration):
        """Генерирует ключ кэша для фрагмента"""
        # Округляем позицию до 1 секунды для группировки близких фрагментов
        rounded_pos = int(position)
        return f"{rounded_pos}_{duration}"
        
    def get_cached_fragment(self, position, duration):
        """Получает фрагмент из кэша"""
        with self.cache_lock:
            key = self.get_cache_key(position, duration)
            return self.cache.get(key)
            
    def cache_fragment(self, position, duration, data):
        """Сохраняет фрагмент в кэш"""
        with self.cache_lock:
            key = self.get_cache_key(position, duration)
            # Ограничиваем размер кэша
            if len(self.cache) > 10:
                # Удаляем самый старый элемент
                oldest_key = next(iter(self.cache))
                del self.cache[oldest_key]
            self.cache[key] = data
            
    def process_fragment_async(self, position, duration):
        """Асинхронная обработка фрагмента в отдельном потоке"""
        if self.is_processing:
            return
            
        self.is_processing = True
        
        def process():
            try:
                if self.mm is None:
                    return
                    
                # Проверяем кэш
                cached = self.get_cached_fragment(position, duration)
                if cached is not None:
                    self.display_queue.put(('cached', cached))
                    return
                
                # Обрабатываем фрагмент
                start_sample = int(position * FS)
                end_sample = min(start_sample + int(duration * FS), self.total)
                
                if end_sample <= start_sample:
                    return
                    
                # Извлекаем данные
                fragment = extract_channel(self.mm, self.total, 0, start_sample, end_sample - start_sample)
                fragment = fragment.astype(np.float64) * ((2*2.4)/(2**24)) * 1e3
                
                # Применяем фильтры
                fragment = highpass_filter(fragment, FS, cutoff=0.5)
                if self.notch_enabled:
                    fragment = notch_filter(fragment, FS)
                
                # Находим пики в фрагменте
                peaks, _, _, _, _ = detect_r_peaks(
                    fragment, FS,
                    sensitivity=self.sensitivity.get(),
                    min_rr_ms=self.min_rr_ms.get(),
                    window_ms=self.window_ms.get(),
                    search_window_ms=self.search_window_ms.get(),
                    notch_enabled=False  # Уже применен
                )
                
                # Корректируем позиции пиков
                peaks += start_sample
                
                # Фильтруем пики в пределах фрагмента
                peaks = peaks[(peaks >= start_sample) & (peaks < end_sample)]
                
                # Конвертируем в относительные позиции для отображения
                display_peaks = [(p - start_sample, fragment[p - start_sample]) for p in peaks]
                
                result = {
                    'fragment': fragment,
                    'peaks': display_peaks,
                    'start_sample': start_sample,
                    'position': position
                }
                
                # Кэшируем результат
                self.cache_fragment(position, duration, result)
                
                # Отправляем в очередь отображения
                self.display_queue.put(('processed', result))
                
            except Exception as e:
                print(f"Ошибка обработки фрагмента: {e}")
            finally:
                self.is_processing = False
        
        # Запускаем в отдельном потоке
        self.executor.submit(process)
    
    def get_current_view_range(self):
        """Получает текущий диапазон просмотра из matplotlib"""
        xlim = self.ax_raw.get_xlim()
        return xlim[0], xlim[1]
        
    def get_current_position(self):
        """Получает текущую позицию просмотра"""
        xlim = self.ax_raw.get_xlim()
        return (xlim[0] + xlim[1]) / 2
    

    
    def update_display(self):
        """Обновление отображения - теперь использует matplotlib навигацию"""
        # Этот метод больше не нужен, так как навигация через matplotlib
        pass

    def draw_raw(self):
        """Отрисовка полной записи с пиками"""
        if not hasattr(self, 'full_ecg') or self.full_ecg is None:
            return
            
        # Сохраняем текущие границы отображения перед очисткой
        current_xlim = self.ax_raw.get_xlim()
        current_ylim = self.ax_raw.get_ylim()
        
        # Проверяем, нужно ли сохранять границы
        # Сохраняем только если это не первая отрисовка и границы были изменены пользователем
        preserve_limits = self._first_draw_complete
        
        # Очищаем график
        self.ax_raw.cla()
        
        # Отрисовываем полную запись
        self.ax_raw.plot(self.t_full, self.full_ecg, color='black', linewidth=0.5, alpha=0.8)
        
        # Отрисовываем пики
        if hasattr(self, 'r_peaks') and len(self.r_peaks) > 0:
            # Проверяем, что индексы пиков не выходят за границы
            valid_peaks = self.r_peaks[(self.r_peaks >= 0) & (self.r_peaks < len(self.full_ecg))]
            
            if len(valid_peaks) > 0:
                peak_times = valid_peaks / FS
                peak_values = self.full_ecg[valid_peaks]
                
                # Разделяем автоматические и ручные пики
                auto_peaks = []
                manual_peaks = []
                selected_peak = None
                
                for i, (t, val) in enumerate(zip(peak_times, peak_values)):
                    # Находим соответствующий индекс в оригинальном массиве пиков
                    orig_idx = np.where(self.r_peaks == valid_peaks[i])[0][0]
                    
                    if orig_idx in self.manual_peaks:
                        manual_peaks.append((t, val))
                    elif orig_idx == self.selected_idx:
                        selected_peak = (t, val)
                    else:
                        auto_peaks.append((t, val))
                
                # Отрисовываем пики
                if auto_peaks:
                    x, y = zip(*auto_peaks)
                    self.ax_raw.plot(x, y, 'ro', markersize=3, label='Автоматические пики')
                
                if manual_peaks:
                    x, y = zip(*manual_peaks)
                    self.ax_raw.plot(x, y, 'go', markersize=4, label='Ручные пики')
                
                if selected_peak:
                    x, y = selected_peak
                    self.ax_raw.plot(x, y, 'bo', markersize=6, label='Выбранный пик')
                
                self.ax_raw.legend(loc='upper right', fontsize=8)
        
        # Отрисовываем выделенные регионы
        for patch in self.region_patches:
            self.ax_raw.add_patch(patch)
        
        # Настройки графика
        self.ax_raw.set_xlabel('Время (сек)')
        self.ax_raw.set_ylabel('Амплитуда (мВ)')
        self.ax_raw.set_title('Полная запись ЭКГ')
        self.ax_raw.grid(True, alpha=0.3)
        
        # Восстанавливаем границы отображения, если они были установлены пользователем
        if preserve_limits:
            self.ax_raw.set_xlim(current_xlim)
            self.ax_raw.set_ylim(current_ylim)
        
        # Отмечаем, что первая отрисовка завершена
        self._first_draw_complete = True
        
        # Обновляем canvas
        self.canvas_raw.draw()

    def draw_selected(self):
        self.ax_sel.cla()
        if self.selected_idx is not None and self.selected_idx < len(self.r_peaks):
            r = self.r_peaks[self.selected_idx]
            half = int(0.3 * FS)
            start = max(r-half, 0)
            end = min(r+half, self.total)
            
            # Извлекаем данные с полной частотой дискретизации
            raw = extract_channel(self.mm, self.total, 0, start, end - start).astype(np.float64)
            seg = raw * ((2*2.4)/(2**24)) * 1e3
            
            # Применяем фильтры
            seg = highpass_filter(seg, FS, cutoff=0.5)
            if self.notch_enabled:
                seg = notch_filter(seg, FS)
            
            tseg = np.linspace(-half/FS, half/FS, len(seg))
            self.ax_sel.plot(tseg, seg, color='green')
        self.ax_sel.set_title('Выбранный комплекс')
        self.ax_sel.set_xlabel('Время (с)')
        self.ax_sel.set_ylabel('Амплитуда (мВ)')
        self.ax_sel.grid(True, alpha=0.3)
        self.canvas_sel.draw()

    def copy_screenshot(self):
        x,y = self.winfo_rootx(), self.winfo_rooty()
        w,h = self.winfo_width(), self.winfo_height()
        img = ImageGrab.grab((x,y,x+w,y+h))
        copy_image_to_clipboard(img)
        messagebox.showinfo('Copied','Снимок в буфере обмена')

if __name__ == '__main__':
    app = ECGApp()
    app.mainloop()




