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

def notch_filter(data, fs, freq=50.0, Q=30.0):
    b, a = iirnotch(freq, Q, fs)
    return filtfilt(b, a, data)

def detect_r_peaks(ecg, fs, sensitivity=1.0, min_rr_ms=200, window_ms=120, search_window_ms=50):
    """
    Алгоритм Томсона (Pan-Tompkins) для детекции R-пиков в ЭКГ
    
    Параметры:
    - sensitivity: чувствительность (0.1-2.0, по умолчанию 1.0)
    - min_rr_ms: минимальное расстояние между пиками в мс (по умолчанию 200)
    - window_ms: размер окна интегрирования в мс (по умолчанию 120)
    - search_window_ms: окно поиска точного пика в мс (по умолчанию 50)
    """
    # 1. Полосовой фильтр (5-15 Гц для QRS)
    from scipy.signal import butter, filtfilt
    
    # Нормализация частоты дискретизации
    nyquist = fs / 2
    low = 5 / nyquist
    high = 15 / nyquist
    b, a = butter(2, [low, high], btype='band')
    filtered = filtfilt(b, a, ecg)
    
    # 2. Дифференцирование
    diff_signal = np.diff(filtered)
    
    # 3. Возведение в квадрат
    squared = diff_signal ** 2
    
    # 4. Интегрирование с окном
    window_size = int(window_ms * fs / 1000)
    integrated = np.convolve(squared, np.ones(window_size), mode='same')
    
    # 5. Адаптивный порог с учетом чувствительности
    # Начальный порог
    threshold_i1 = (0.6 / sensitivity) * np.max(integrated)
    threshold_i2 = 0.6 * threshold_i1
    
    # Списки для хранения пиков
    r_peaks = []
    noise_peaks = []
    signal_peaks = []
    
    # Минимальное расстояние между пиками
    min_rr = int(min_rr_ms * fs / 1000)
    
    # Поиск пиков
    for i in range(1, len(integrated) - 1):
        if integrated[i] > integrated[i-1] and integrated[i] > integrated[i+1]:
            if integrated[i] > threshold_i1:
                # Проверяем минимальное расстояние
                if not r_peaks or (i - r_peaks[-1]) > min_rr:
                    r_peaks.append(i)
                    signal_peaks.append(i)
            elif integrated[i] > threshold_i2:
                noise_peaks.append(i)
    
    # 6. Адаптивное обновление порогов
    if len(signal_peaks) > 8:
        # Обновляем пороги каждые 8 пиков
        signal_level = np.mean([integrated[peak] for peak in signal_peaks[-8:]])
        noise_level = np.mean([integrated[peak] for peak in noise_peaks[-8:]]) if noise_peaks else signal_level * 0.5
        
        threshold_i1 = (noise_level + 0.25 * (signal_level - noise_level)) / sensitivity
        threshold_i2 = 0.5 * threshold_i1
    
    # 7. Поиск точных позиций R-пиков в оригинальном сигнале
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

def copy_image_to_clipboard(img):
    output = io.BytesIO()
    img.convert("RGB").save(output, "BMP")
    data = output.getvalue()[14:]
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32con.CF_DIB, data)
    win32clipboard.CloseClipboard()

def plot_tompson_stages(ax, ecg, peaks, filtered, diff_signal, squared, integrated, fs):
    """
    Отображение этапов алгоритма Томсона
    """
    ax.clear()
    
    # Показываем только первые 5 секунд для наглядности
    t_max = 5.0
    n_samples = int(t_max * fs)
    t = np.arange(min(n_samples, len(ecg))) / fs
    
    # Нормализация сигналов для отображения
    ecg_norm = ecg[:len(t)] / np.max(np.abs(ecg[:len(t)]))
    filtered_norm = filtered[:len(t)] / np.max(np.abs(filtered[:len(t)])) if len(filtered) >= len(t) else filtered / np.max(np.abs(filtered))
    diff_norm = diff_signal[:len(t)-1] / np.max(np.abs(diff_signal[:len(t)-1])) if len(diff_signal) >= len(t)-1 else diff_signal / np.max(np.abs(diff_signal))
    squared_norm = squared[:len(t)-1] / np.max(squared[:len(t)-1]) if len(squared) >= len(t)-1 else squared / np.max(squared)
    integrated_norm = integrated[:len(t)] / np.max(integrated[:len(t)]) if len(integrated) >= len(t) else integrated / np.max(integrated)
    
    # Смещение для отображения
    offset = 2.5
    ax.plot(t, ecg_norm + offset*4, 'b-', linewidth=0.8, label='Исходный ЭКГ')
    ax.plot(t, filtered_norm + offset*3, 'g-', linewidth=0.8, label='Полосовой фильтр')
    ax.plot(t[:-1], diff_norm + offset*2, 'r-', linewidth=0.8, label='Дифференцирование')
    ax.plot(t[:-1], squared_norm + offset*1, 'm-', linewidth=0.8, label='Квадрат')
    ax.plot(t, integrated_norm, 'k-', linewidth=0.8, label='Интегрирование')
    
    # Отмечаем найденные пики
    peaks_in_window = peaks[peaks < n_samples]
    if len(peaks_in_window) > 0:
        ax.plot(peaks_in_window/fs, integrated_norm[peaks_in_window], 'ro', markersize=4, label='R-пики')
    
    ax.set_title('Этапы алгоритма Томсона')
    ax.set_ylabel('Амплитуда (норм.)')
    ax.set_xlabel('Время (с)')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)

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
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Открыть запись...", command=self.open_file)
        filemenu.add_separator()
        filemenu.add_command(label="Сохранить параметры...", command=self.save_params)
        filemenu.add_command(label="Загрузить параметры...", command=self.load_params)
        filemenu.add_separator()
        filemenu.add_command(label="Выход", command=lambda: on_close(self))
        menubar.add_cascade(label="Файл", menu=filemenu)
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
        self.select_beats = False
        self.select_complex_mode = False
        self.selected_idx = None
        
        # --- Параметры алгоритма Томсона ---
        self.sensitivity = tk.DoubleVar(value=1.0)
        self.min_rr_ms = tk.IntVar(value=200)
        self.window_ms = tk.IntVar(value=120)
        self.search_window_ms = tk.IntVar(value=50)

        # --- Разметка grid ---
        self.rowconfigure(0, weight=3)
        self.rowconfigure(1, weight=2)
        self.rowconfigure(2, weight=2)
        self.rowconfigure(3, weight=2)
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

        # --- Спектр ---
        self.fig_sp, self.ax_sp = plt.subplots(figsize=(4,2))
        plt.tight_layout()
        self.canvas_sp = FigureCanvasTkAgg(self.fig_sp, master=self)
        self.canvas_sp.get_tk_widget().grid(row=1, column=1, sticky='nsew', padx=5, pady=5)

        # --- Этапы алгоритма Томсона ---
        self.fig_tom, self.ax_tom = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_tom = FigureCanvasTkAgg(self.fig_tom, master=self)
        self.canvas_tom.get_tk_widget().grid(row=2, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

        # --- Выбранный комплекс ---
        self.fig_sel, self.ax_sel = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_sel = FigureCanvasTkAgg(self.fig_sel, master=self)
        self.canvas_sel.get_tk_widget().grid(row=3, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

        # --- Панель управления алгоритмом Томсона ---
        tompson_frame = ttk.LabelFrame(self, text="Параметры алгоритма Томсона")
        tompson_frame.grid(row=4, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        
        # Первая строка параметров
        param_row1 = ttk.Frame(tompson_frame)
        param_row1.pack(fill='x', padx=5, pady=2)
        
        ttk.Label(param_row1, text="Чувствительность:").pack(side='left', padx=5)
        sensitivity_scale = ttk.Scale(param_row1, from_=0.1, to=2.0, variable=self.sensitivity, 
                                     orient='horizontal', length=150, command=self.on_param_change)
        sensitivity_scale.pack(side='left', padx=5)
        self.sensitivity_label = ttk.Label(param_row1, text="1.0")
        self.sensitivity_label.pack(side='left', padx=5)
        
        ttk.Label(param_row1, text="Мин. RR (мс):").pack(side='left', padx=5)
        min_rr_scale = ttk.Scale(param_row1, from_=100, to=500, variable=self.min_rr_ms, 
                                orient='horizontal', length=150, command=self.on_param_change)
        min_rr_scale.pack(side='left', padx=5)
        self.min_rr_label = ttk.Label(param_row1, text="200")
        self.min_rr_label.pack(side='left', padx=5)
        
        # Вторая строка параметров
        param_row2 = ttk.Frame(tompson_frame)
        param_row2.pack(fill='x', padx=5, pady=2)
        
        ttk.Label(param_row2, text="Окно интегрирования (мс):").pack(side='left', padx=5)
        window_scale = ttk.Scale(param_row2, from_=50, to=200, variable=self.window_ms, 
                                orient='horizontal', length=150, command=self.on_param_change)
        window_scale.pack(side='left', padx=5)
        self.window_label = ttk.Label(param_row2, text="120")
        self.window_label.pack(side='left', padx=5)
        
        ttk.Label(param_row2, text="Окно поиска (мс):").pack(side='left', padx=5)
        search_scale = ttk.Scale(param_row2, from_=20, to=100, variable=self.search_window_ms, 
                                orient='horizontal', length=150, command=self.on_param_change)
        search_scale.pack(side='left', padx=5)
        self.search_label = ttk.Label(param_row2, text="50")
        self.search_label.pack(side='left', padx=5)
        
        # --- Панель кнопок ---
        ctrl = ttk.Frame(self)
        ctrl.grid(row=5, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        ttk.Button(ctrl, text='Notch 50 Hz', command=self.toggle_notch).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Выделить регион', command=self.activate_region_selector).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Отметить комплекс', command=self.toggle_peak_selector).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Выбрать комплекс', command=self.toggle_select_complex).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Обновить', command=self.compute).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Сброс параметров', command=self.reset_tompson_params).pack(side='left', padx=5)
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
            self.span.set_active(True)

    def toggle_peak_selector(self):
        self.select_beats = not self.select_beats

    def toggle_select_complex(self):
        self.select_complex_mode = not self.select_complex_mode

    def on_param_change(self, event=None):
        """Обработчик изменения параметров алгоритма Томсона"""
        # Обновляем метки
        self.sensitivity_label.config(text=f"{self.sensitivity.get():.1f}")
        self.min_rr_label.config(text=f"{self.min_rr_ms.get()}")
        self.window_label.config(text=f"{self.window_ms.get()}")
        self.search_label.config(text=f"{self.search_window_ms.get()}")
        
        # Автоматически пересчитываем при изменении параметров
        if self.mm is not None:
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

    def clear_markup(self):
        # удалить патчи
        for p in self.region_patches:
            try: p.remove()
            except: pass
        self.region_patches = []
        self.bad_regions = []
        self.bad_beats = set()
        self.manual_beats = set()
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
        # правка пиков
        if self.select_beats and ev.button==3:
            t = ev.xdata
            idx = np.argmin(np.abs(self.r_peaks/FS - t))
            if abs(self.r_peaks[idx]/FS - t) < 0.05:
                if idx in self.bad_beats: self.bad_beats.remove(idx)
                else: self.bad_beats.add(idx)
                self.draw_raw()
        # выбор комплекса
        if self.select_complex_mode and ev.button==1:
            idx = np.argmin(np.abs(self.r_peaks/FS - ev.xdata))
            self.selected_idx = idx
            self.draw_selected()

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
        if self.notch_enabled:
            ecg = notch_filter(ecg, FS)
        self.ecg = ecg

        self.r_peaks, filtered, diff_signal, squared, integrated = detect_r_peaks(
            ecg, FS, 
            sensitivity=self.sensitivity.get(),
            min_rr_ms=self.min_rr_ms.get(),
            window_ms=self.window_ms.get(),
            search_window_ms=self.search_window_ms.get()
        )
        # HR по годным пикам
        good = [r for i,r in enumerate(self.r_peaks)
                if i not in self.bad_beats and not any(a<=r/FS<=b for a,b in self.bad_regions)]
        # Обновляем информацию о пиках
        self.peaks_label.config(text=f"Пики: {len(self.r_peaks)}")
        
        if len(good)>1:
            rr = np.diff(np.array(good)/FS)
            hr = (60/rr).mean()
            self.hr_label.config(text=f"ЧСС: {hr:.0f} bpm")
        else:
            self.hr_label.config(text="ЧСС: -- bpm")

        self.draw_raw()
        
        allb = np.unique(np.concatenate([good, list(self.manual_beats)])).astype(int)
        t_avg, avg = average_complex(ecg, allb, FS)

        self.ax_avg.cla()
        if avg is not None:
            self.ax_avg.plot(t_avg, avg, linewidth=1)
        self.ax_avg.set_title('Усреднённый комплекс')
        self.canvas_avg.draw()

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

        # Отображение этапов алгоритма Томсона
        plot_tompson_stages(self.ax_tom, ecg, self.r_peaks, filtered, diff_signal, squared, integrated, FS)
        self.canvas_tom.draw()

        self.draw_selected()

    def draw_raw(self):
        self.ax_raw.cla()
        t = np.arange(len(self.ecg)) / FS
        self.ax_raw.plot(t, self.ecg, color='black', linewidth=0.8)
        self.ax_raw.plot(self.r_peaks/FS, self.ecg[self.r_peaks], 'ro', markersize=3)
        for patch in self.region_patches:
            self.ax_raw.add_patch(patch)
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




