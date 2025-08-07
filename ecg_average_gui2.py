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

# === Параметры сигнала ===
FS = 2000.0
BYTES_ECG = 3
BYTES_AUX = 2 * 4
N_CHANNELS = 8
SAMPLE_SIZE = BYTES_ECG * N_CHANNELS + BYTES_AUX
START_DELAY = 0.5
READ_COUNT = int(100 * FS)

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

def create_ecg_markers(t, signal):
    """
    Создание фиксированных маркеров ЭКГ для ручного позиционирования
    """
    markers = {}
    
    # Находим примерный центр комплекса (максимум сигнала)
    center_idx = np.argmax(signal)
    center_time = t[center_idx]
    
    # Создаем маркеры с примерными позициями
    # Начало P (примерно 200 мс перед центром)
    p_start_idx = max(0, np.where(t >= center_time - 0.2)[0][0] if len(np.where(t >= center_time - 0.2)[0]) > 0 else 0)
    markers['начало P'] = (t[p_start_idx], signal[p_start_idx])
    
    # Конец P (примерно 80 мс перед центром)
    p_end_idx = max(0, np.where(t >= center_time - 0.08)[0][0] if len(np.where(t >= center_time - 0.08)[0]) > 0 else 0)
    markers['конец P'] = (t[p_end_idx], signal[p_end_idx])
    
    # Q (примерно 40 мс перед центром)
    q_idx = max(0, np.where(t >= center_time - 0.04)[0][0] if len(np.where(t >= center_time - 0.04)[0]) > 0 else 0)
    markers['Q'] = (t[q_idx], signal[q_idx])
    
    # R (центр)
    markers['R'] = (t[center_idx], signal[center_idx])
    
    # S (примерно 40 мс после центра)
    s_idx = min(len(signal)-1, np.where(t >= center_time + 0.04)[0][0] if len(np.where(t >= center_time + 0.04)[0]) > 0 else len(signal)-1)
    markers['S'] = (t[s_idx], signal[s_idx])
    
    # Начало T (примерно 80 мс после центра)
    t_start_idx = min(len(signal)-1, np.where(t >= center_time + 0.08)[0][0] if len(np.where(t >= center_time + 0.08)[0]) > 0 else len(signal)-1)
    markers['начало T'] = (t[t_start_idx], signal[t_start_idx])
    
    # Конец T (примерно 200 мс после центра)
    t_end_idx = min(len(signal)-1, np.where(t >= center_time + 0.2)[0][0] if len(np.where(t >= center_time + 0.2)[0]) > 0 else len(signal)-1)
    markers['конец T'] = (t[t_end_idx], signal[t_end_idx])
    
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
        
        # --- Параметры алгоритма Томсона ---
        self.sensitivity = tk.DoubleVar(value=1.0)
        self.min_rr_ms = tk.IntVar(value=200)
        self.window_ms = tk.IntVar(value=120)
        self.search_window_ms = tk.IntVar(value=50)
        self.auto_recompute = tk.BooleanVar(value=False)  # Контроль автоматического пересчета

        # --- Разметка grid ---
        self.rowconfigure(0, weight=3)
        self.rowconfigure(1, weight=2)
        self.rowconfigure(2, weight=2)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(4, weight=1)
        self.rowconfigure(5, weight=0)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        # --- Raw ECG и HR ---
        self.fig_raw, self.ax_raw = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_raw = FigureCanvasTkAgg(self.fig_raw, master=self)
        self.canvas_raw.get_tk_widget().grid(row=0, column=0, columnspan=2, sticky='nsew')
        self.hr_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="ЧСС: -- bpm", background='#fff')
        self.hr_label.place(relx=0.01, rely=0.01)
        
        self.peaks_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="Пики: 0", background='#fff')
        self.peaks_label.place(relx=0.01, rely=0.05)
        
        self.progress_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="", background='#fff')
        self.progress_label.place(relx=0.01, rely=0.09)

        # SpanSelector и события мыши
        self.span = SpanSelector(self.ax_raw, self.onselect_region, 'horizontal', useblit=True,
                                 props=dict(alpha=0.3, facecolor='red'), button=1, minspan=0.01)
        self.span.set_active(False)
        self.fig_raw.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig_raw.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig_raw.canvas.mpl_connect('button_release_event', self.on_release)

        # --- Усреднённый комплекс ---
        self.fig_avg, self.ax_avg = plt.subplots(figsize=(4,2))
        plt.tight_layout()
        self.canvas_avg = FigureCanvasTkAgg(self.fig_avg, master=self)
        self.canvas_avg.get_tk_widget().grid(row=1, column=0, sticky='nsew', padx=5, pady=5)
        
        # Подключаем обработчик кликов для усредненного комплекса
        self.fig_avg.canvas.mpl_connect('button_press_event', self.on_avg_click)

        # --- Спектр ---
        self.fig_sp, self.ax_sp = plt.subplots(figsize=(4,2))
        plt.tight_layout()
        self.canvas_sp = FigureCanvasTkAgg(self.fig_sp, master=self)
        self.canvas_sp.get_tk_widget().grid(row=1, column=1, sticky='nsew', padx=5, pady=5)



        # --- Выбранный комплекс ---
        self.fig_sel, self.ax_sel = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_sel = FigureCanvasTkAgg(self.fig_sel, master=self)
        self.canvas_sel.get_tk_widget().grid(row=2, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

        # --- Текстовое поле для параметров ---
        text_frame = ttk.LabelFrame(self, text="Параметры и измерения")
        text_frame.grid(row=3, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        
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
        ctrl.grid(row=4, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        ttk.Button(ctrl, text='Notch 50 Hz', command=self.toggle_notch).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Выделить регион', command=self.activate_region_selector).pack(side='left', padx=5)
        
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
        except Exception as e:
            messagebox.showerror('Ошибка', str(e))

    def toggle_notch(self):
        self.notch_enabled = not self.notch_enabled
        self.compute()

    def activate_region_selector(self):
        if self.mm is not None:
            # Отключаем другие режимы
            self.manual_peak_mode = False
            self.select_complex_mode = False
            
            # Обновляем состояние кнопок
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            
            # Активируем SpanSelector
            self.span.set_active(True)

    def toggle_peak_selector(self):
        """Переключение режима ручного добавления/удаления пиков"""
        self.manual_peak_mode = not self.manual_peak_mode
        self.select_complex_mode = False  # отключаем другой режим
        
        # Обновляем состояние кнопок
        if self.manual_peak_mode:
            self.manual_peak_button.config(bg='lightgreen', relief='sunken')
            self.select_complex_button.config(bg='lightgray', relief='raised')
        else:
            self.manual_peak_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector при активном режиме ручной разметки
        self.span.set_active(False)

    def toggle_select_complex(self):
        """Переключение режима выбора комплекса для отображения"""
        self.select_complex_mode = not self.select_complex_mode
        self.manual_peak_mode = False  # отключаем другой режим
        
        # Обновляем состояние кнопок
        if self.select_complex_mode:
            self.select_complex_button.config(bg='lightblue', relief='sunken')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
        else:
            self.select_complex_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector при активном режиме выбора комплекса
        self.span.set_active(False)

    def toggle_measurement_mode(self):
        """Переключение режима измерений на усредненном комплексе"""
        self.measurement_mode = not self.measurement_mode
        self.manual_peak_mode = False
        self.select_complex_mode = False
        
        # Обновляем состояние кнопок
        if self.measurement_mode:
            self.measurement_button.config(bg='lightyellow', relief='sunken')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            messagebox.showinfo("Режим измерений", 
                              "Включен режим измерений.\n"
                              "Кликните на усредненном комплексе для добавления маркеров.\n"
                              "Правый клик - удалить маркер")
        else:
            self.measurement_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector
        self.span.set_active(False)

    def toggle_marker_edit(self):
        """Переключение режима редактирования маркеров ЭКГ"""
        self.marker_edit_mode = not self.marker_edit_mode
        self.measurement_mode = False
        self.manual_peak_mode = False
        self.select_complex_mode = False
        
        # Обновляем состояние кнопок
        if self.marker_edit_mode:
            self.marker_edit_button.config(bg='lightcoral', relief='sunken')
            self.measurement_button.config(bg='lightgray', relief='raised')
            self.manual_peak_button.config(bg='lightgray', relief='raised')
            self.select_complex_button.config(bg='lightgray', relief='raised')
            messagebox.showinfo("Редактирование маркеров", 
                              "Включен режим редактирования маркеров ЭКГ.\n"
                              "Кликните на маркер для его перемещения.\n"
                              "Правый клик - удалить маркер")
        else:
            self.marker_edit_button.config(bg='lightgray', relief='raised')
        
        # Отключаем SpanSelector
        self.span.set_active(False)

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
        if self.ecg is None or len(self.r_peaks) == 0:
            messagebox.showwarning("Предупреждение", "Нет данных для анализа")
            return
        
        # Вычисляем параметры
        good_peaks = [r for r in self.r_peaks if not any(a<=r/FS<=b for a,b in self.bad_regions)]
        
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
                marker_order = ['начало P', 'конец P', 'Q', 'R', 'S', 'начало T', 'конец T']
                for marker_name in marker_order:
                    if marker_name in all_markers:
                        x, y = all_markers[marker_name]
                        status = "(ручная)" if marker_name in self.manual_markers else "(авто)"
                        text_content += f"• {marker_name} {status}: {x:.1f} мс, {y:.2f} мВ\n"
                
                # Вычисляем интервалы
                if 'конец P' in all_markers and 'Q' in all_markers:
                    pq_interval = (all_markers['Q'][0] - all_markers['конец P'][0]) * 1000  # в мс
                    text_content += f"• PQ интервал: {pq_interval:.1f} мс\n"
                
                if 'Q' in all_markers and 'S' in all_markers:
                    qrs_duration = (all_markers['S'][0] - all_markers['Q'][0]) * 1000  # в мс
                    text_content += f"• QRS длительность: {qrs_duration:.1f} мс\n"
                
                if 'S' in all_markers and 'начало T' in all_markers:
                    st_interval = (all_markers['начало T'][0] - all_markers['S'][0]) * 1000  # в мс
                    text_content += f"• ST интервал: {st_interval:.1f} мс\n"
                
                if 'начало P' in all_markers and 'конец P' in all_markers:
                    p_duration = (all_markers['конец P'][0] - all_markers['начало P'][0]) * 1000  # в мс
                    text_content += f"• Длительность P: {p_duration:.1f} мс\n"
                
                if 'начало T' in all_markers and 'конец T' in all_markers:
                    t_duration = (all_markers['конец T'][0] - all_markers['начало T'][0]) * 1000  # в мс
                    text_content += f"• Длительность T: {t_duration:.1f} мс\n"
            
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
                # Добавляем точку измерения
                self.measurement_points.append((event.xdata, event.ydata, f"Точка {len(self.measurement_points)+1}"))
                
                # Отображаем маркер
                self.ax_avg.plot(event.xdata, event.ydata, 'ro', markersize=6)
                self.ax_avg.annotate(f"Точка {len(self.measurement_points)}", 
                                   (event.xdata, event.ydata), 
                                   xytext=(5, 5), textcoords='offset points',
                                   fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
                
                self.canvas_avg.draw()
                
            elif event.button == 3:  # Правый клик - удалить ближайший маркер
                if self.measurement_points:
                    # Находим ближайшую точку
                    distances = [(abs(event.xdata - x), i) for i, (x, y, _) in enumerate(self.measurement_points)]
                    distances.sort()
                    nearest_idx = distances[0][1]
                    
                    # Удаляем точку
                    del self.measurement_points[nearest_idx]
                    
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
                    
                    # Перемещаем маркер
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
        
        if avg is not None:
            self.ax_avg.plot(t_avg, avg, linewidth=1, color='blue')
            
            # Отображаем маркеры измерений
            for i, (x, y, label) in enumerate(self.measurement_points, 1):
                self.ax_avg.plot(x, y, 'ro', markersize=6)
                self.ax_avg.annotate(f"Точка {i}", 
                                   (x, y), 
                                   xytext=(5, 5), textcoords='offset points',
                                   fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
            
            # Отображаем маркеры ЭКГ
            self.draw_ecg_markers(t_avg, avg)
        
        self.ax_avg.set_title('Усреднённый комплекс')
        self.ax_avg.set_xlabel('Время (мс)')
        self.ax_avg.set_ylabel('Амплитуда (мВ)')
        self.ax_avg.grid(True, alpha=0.3)
        self.canvas_avg.draw()

    def draw_ecg_markers(self, t_avg, avg):
        """Отрисовка маркеров ЭКГ"""
        # Цвета для разных маркеров
        colors = {
            'начало P': 'lightgreen', 'конец P': 'green',
            'Q': 'orange', 'R': 'red', 'S': 'purple',
            'начало T': 'lightblue', 'конец T': 'blue'
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
                                   fontsize=8, fontweight='bold', color=color,
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
        self.compute()

    def onselect_region(self, x0, x1):
        start, end = min(x0,x1), max(x0,x1)
        self.bad_regions.append((start, end))
        patch = Rectangle((start, min(self.ecg)), end-start, max(self.ecg)-min(self.ecg),
                          color='red', alpha=0.3)
        self.ax_raw.add_patch(patch)
        self.region_patches.append(patch)
        self.canvas_raw.draw()

    def on_click(self, ev):
        if ev.inaxes != self.ax_raw:
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
            t = ev.xdata
            sample_idx = int(t * FS)
            
            if ev.button == 3:  # Правый клик - добавить пик
                # Ищем локальный максимум в окрестности клика
                search_window = int(0.1 * FS)  # 100 мс окно поиска
                start = max(0, sample_idx - search_window)
                end = min(len(self.ecg), sample_idx + search_window)
                
                # Находим локальный максимум
                local_max_idx = start + np.argmax(self.ecg[start:end])
                
                # Добавляем пик в список
                self.r_peaks = np.append(self.r_peaks, local_max_idx)
                self.r_peaks = np.sort(self.r_peaks)  # сортируем
                # Находим новый индекс добавленного пика
                new_idx = np.where(self.r_peaks == local_max_idx)[0][0]
                self.manual_peaks.add(new_idx)  # запоминаем как ручной
                
                self.draw_raw()
                
            elif ev.button == 1:  # Левый клик - удалить пик
                # Ищем ближайший пик
                if len(self.r_peaks) > 0:
                    distances = np.abs(self.r_peaks/FS - t)
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
                        
                        self.draw_raw()
        
        # Выбор комплекса для отображения
        elif self.select_complex_mode and ev.button == 1:
            if len(self.r_peaks) > 0:
                # Ищем ближайший пик
                distances = np.abs(self.r_peaks/FS - ev.xdata)
                nearest_idx = np.argmin(distances)
                
                if distances[nearest_idx] < 0.1:  # 100 мс допуск
                    self.selected_idx = nearest_idx
                    self.draw_selected()
                    self.draw_raw()  # Обновляем отображение для показа выбранного пика

    def on_motion(self, ev):
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
        self.resizing = None

    def compute(self):
        if self.mm is None:
            return
        raw = extract_channel(self.mm, self.total, 0, int(START_DELAY*FS), READ_COUNT).astype(np.float64)
        ecg = raw * ((2*2.4)/(2**24)) * 1e3
        
        # Применяем фильтр высоких частот для устранения плавающей изолинии
        ecg = highpass_filter(ecg, FS, cutoff=0.5)
        
        if self.notch_enabled:
            ecg = notch_filter(ecg, FS)
        self.ecg = ecg

        self.r_peaks, filtered, diff_signal, squared, integrated = detect_r_peaks(
            ecg, FS, 
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

        self.draw_raw()
        
        t_avg, avg = average_complex(ecg, good, FS)
        
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

    def draw_raw(self):
        self.ax_raw.cla()
        t = np.arange(len(self.ecg)) / FS
        self.ax_raw.plot(t, self.ecg, color='black', linewidth=0.8)
        
        # Отображаем автоматически найденные пики красным
        auto_peaks = [i for i, peak in enumerate(self.r_peaks) if i not in self.manual_peaks]
        if auto_peaks:
            auto_peak_positions = self.r_peaks[auto_peaks]
            self.ax_raw.plot(auto_peak_positions/FS, self.ecg[auto_peak_positions], 'ro', markersize=3, label='Автоматические пики')
        
        # Отображаем ручные пики зеленым
        if self.manual_peaks:
            manual_peak_positions = self.r_peaks[list(self.manual_peaks)]
            self.ax_raw.plot(manual_peak_positions/FS, self.ecg[manual_peak_positions], 'go', markersize=4, label='Ручные пики')
        
        # Отображаем выбранный пик синим
        if self.selected_idx is not None and self.selected_idx < len(self.r_peaks):
            selected_peak = self.r_peaks[self.selected_idx]
            self.ax_raw.plot(selected_peak/FS, self.ecg[selected_peak], 'bo', markersize=6, label='Выбранный пик')
        
        for patch in self.region_patches:
            self.ax_raw.add_patch(patch)
        
        if len(self.r_peaks) > 0:
            self.ax_raw.legend(loc='upper right', fontsize=8)
        
        self.canvas_raw.draw()

    def draw_selected(self):
        self.ax_sel.cla()
        if self.selected_idx is not None:
            r = self.r_peaks[self.selected_idx]
            half = int(0.3 * FS)
            seg = self.ecg[max(r-half,0):r+half]
            tseg = np.linspace(-half/FS, half/FS, len(seg))
            self.ax_sel.plot(tseg, seg, color='green')
        self.ax_sel.set_title('Выбранный комплекс')
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




