# какая версия пуша?
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirnotch, filtfilt, find_peaks, butter, sosfiltfilt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
from matplotlib.widgets import SpanSelector
from PIL import ImageGrab
import io
import win32clipboard
import win32con
import json
import os
from datetime import datetime
import pandas as pd

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

def butter_bandpass_filter(data, fs, lowcut=0.5, highcut=100.0, order=4):
    """Полосовой фильтр Баттерворта"""
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], btype='band', output='sos')
    return sosfiltfilt(sos, data)

def butter_lowpass_filter(data, fs, cutoff=100.0, order=4):
    """Низкочастотный фильтр Баттерворта"""
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    sos = butter(order, normal_cutoff, btype='low', output='sos')
    return sosfiltfilt(sos, data)

def butter_highpass_filter(data, fs, cutoff=0.5, order=4):
    """Высокочастотный фильтр Баттерворта"""
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    sos = butter(order, normal_cutoff, btype='high', output='sos')
    return sosfiltfilt(sos, data)

def detect_r_peaks(ecg, fs):
    height = np.mean(ecg) + 0.15 * (ecg.max() - np.mean(ecg))
    min_dist = int(0.2 * fs)
    peaks, _ = find_peaks(ecg, distance=min_dist, height=height)
    return peaks

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

def calculate_rr_statistics(r_peaks, fs):
    """Вычисление статистики RR интервалов"""
    if len(r_peaks) < 2:
        return {}
    
    rr_intervals = np.diff(r_peaks) / fs * 1000  # в миллисекундах
    hr_intervals = 60000 / rr_intervals  # ЧСС в bpm
    
    stats = {
        'mean_rr': np.mean(rr_intervals),
        'std_rr': np.std(rr_intervals),
        'mean_hr': np.mean(hr_intervals),
        'std_hr': np.std(hr_intervals),
        'min_rr': np.min(rr_intervals),
        'max_rr': np.max(rr_intervals),
        'min_hr': np.min(hr_intervals),
        'max_hr': np.max(hr_intervals),
        'rmssd': np.sqrt(np.mean(np.diff(rr_intervals)**2)),  # RMSSD
        'nn50': np.sum(np.abs(np.diff(rr_intervals)) > 50),  # NN50
        'pnn50': np.sum(np.abs(np.diff(rr_intervals)) > 50) / len(rr_intervals) * 100,  # pNN50
        'total_beats': len(r_peaks)
    }
    return stats

def copy_image_to_clipboard(img):
    output = io.BytesIO()
    img.convert("RGB").save(output, "BMP")
    data = output.getvalue()[14:]
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardData(win32con.CF_DIB, data)
    win32clipboard.CloseClipboard()

def on_close(app):
    app.save_settings()
    plt.close('all')
    app.destroy()
    sys.exit(0)

# === Дополнительные окна ===
class SettingsWindow(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Настройки анализа")
        self.geometry("400x500")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        
        # Параметры фильтров
        ttk.Label(self, text="Параметры фильтров", font=('Arial', 12, 'bold')).pack(pady=10)
        
        # Notch фильтр
        notch_frame = ttk.LabelFrame(self, text="Notch фильтр (50 Hz)")
        notch_frame.pack(fill='x', padx=10, pady=5)
        self.notch_freq = tk.DoubleVar(value=50.0)
        self.notch_q = tk.DoubleVar(value=30.0)
        ttk.Label(notch_frame, text="Частота (Гц):").grid(row=0, column=0, padx=5, pady=2)
        ttk.Entry(notch_frame, textvariable=self.notch_freq, width=10).grid(row=0, column=1, padx=5, pady=2)
        ttk.Label(notch_frame, text="Добротность:").grid(row=1, column=0, padx=5, pady=2)
        ttk.Entry(notch_frame, textvariable=self.notch_q, width=10).grid(row=1, column=1, padx=5, pady=2)
        
        # Полосовой фильтр
        band_frame = ttk.LabelFrame(self, text="Полосовой фильтр")
        band_frame.pack(fill='x', padx=10, pady=5)
        self.band_low = tk.DoubleVar(value=0.5)
        self.band_high = tk.DoubleVar(value=100.0)
        self.band_order = tk.IntVar(value=4)
        ttk.Label(band_frame, text="Нижняя частота (Гц):").grid(row=0, column=0, padx=5, pady=2)
        ttk.Entry(band_frame, textvariable=self.band_low, width=10).grid(row=0, column=1, padx=5, pady=2)
        ttk.Label(band_frame, text="Верхняя частота (Гц):").grid(row=1, column=0, padx=5, pady=2)
        ttk.Entry(band_frame, textvariable=self.band_high, width=10).grid(row=1, column=1, padx=5, pady=2)
        ttk.Label(band_frame, text="Порядок фильтра:").grid(row=2, column=0, padx=5, pady=2)
        ttk.Entry(band_frame, textvariable=self.band_order, width=10).grid(row=2, column=1, padx=5, pady=2)
        
        # Параметры детекции
        detect_frame = ttk.LabelFrame(self, text="Детекция R-пиков")
        detect_frame.pack(fill='x', padx=10, pady=5)
        self.peak_height_factor = tk.DoubleVar(value=0.15)
        self.min_peak_distance = tk.DoubleVar(value=0.2)
        ttk.Label(detect_frame, text="Фактор высоты пика:").grid(row=0, column=0, padx=5, pady=2)
        ttk.Entry(detect_frame, textvariable=self.peak_height_factor, width=10).grid(row=0, column=1, padx=5, pady=2)
        ttk.Label(detect_frame, text="Мин. расстояние (с):").grid(row=1, column=0, padx=5, pady=2)
        ttk.Entry(detect_frame, textvariable=self.min_peak_distance, width=10).grid(row=1, column=1, padx=5, pady=2)
        
        # Параметры усреднения
        avg_frame = ttk.LabelFrame(self, text="Усреднение комплекса")
        avg_frame.pack(fill='x', padx=10, pady=5)
        self.avg_window = tk.IntVar(value=700)
        ttk.Label(avg_frame, text="Окно усреднения (мс):").grid(row=0, column=0, padx=5, pady=2)
        ttk.Entry(avg_frame, textvariable=self.avg_window, width=10).grid(row=0, column=1, padx=5, pady=2)
        
        # Кнопки
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', padx=10, pady=10)
        ttk.Button(btn_frame, text="Применить", command=self.apply_settings).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Отмена", command=self.destroy).pack(side='right', padx=5)
    
    def apply_settings(self):
        # Здесь можно добавить валидацию параметров
        self.destroy()

class StatisticsWindow(tk.Toplevel):
    def __init__(self, parent, rr_stats):
        super().__init__(parent)
        self.title("Статистика RR интервалов")
        self.geometry("500x400")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        
        # Создаем текстовое поле с прокруткой
        text_frame = ttk.Frame(self)
        text_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        self.text_widget = tk.Text(text_frame, wrap='word', font=('Courier', 10))
        scrollbar = ttk.Scrollbar(text_frame, orient='vertical', command=self.text_widget.yview)
        self.text_widget.configure(yscrollcommand=scrollbar.set)
        
        self.text_widget.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        # Заполняем статистикой
        self.display_statistics(rr_stats)
        
        # Кнопки
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', padx=10, pady=10)
        ttk.Button(btn_frame, text="Экспорт в CSV", command=lambda: self.export_to_csv(rr_stats)).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Закрыть", command=self.destroy).pack(side='right', padx=5)
    
    def display_statistics(self, stats):
        if not stats:
            self.text_widget.insert('end', "Недостаточно данных для анализа\n")
            return
        
        self.text_widget.insert('end', "СТАТИСТИКА RR ИНТЕРВАЛОВ\n")
        self.text_widget.insert('end', "=" * 40 + "\n\n")
        
        self.text_widget.insert('end', f"Общее количество ударов: {stats['total_beats']}\n\n")
        
        self.text_widget.insert('end', "RR ИНТЕРВАЛЫ (мс):\n")
        self.text_widget.insert('end', f"  Среднее значение: {stats['mean_rr']:.1f}\n")
        self.text_widget.insert('end', f"  Стандартное отклонение: {stats['std_rr']:.1f}\n")
        self.text_widget.insert('end', f"  Минимальное: {stats['min_rr']:.1f}\n")
        self.text_widget.insert('end', f"  Максимальное: {stats['max_rr']:.1f}\n\n")
        
        self.text_widget.insert('end', "ЧАСТОТА СЕРДЦЕБИЕНИЯ (уд/мин):\n")
        self.text_widget.insert('end', f"  Среднее значение: {stats['mean_hr']:.1f}\n")
        self.text_widget.insert('end', f"  Стандартное отклонение: {stats['std_hr']:.1f}\n")
        self.text_widget.insert('end', f"  Минимальное: {stats['min_hr']:.1f}\n")
        self.text_widget.insert('end', f"  Максимальное: {stats['max_hr']:.1f}\n\n")
        
        self.text_widget.insert('end', "ПОКАЗАТЕЛИ ВАРИАБЕЛЬНОСТИ:\n")
        self.text_widget.insert('end', f"  RMSSD: {stats['rmssd']:.1f} мс\n")
        self.text_widget.insert('end', f"  NN50: {stats['nn50']}\n")
        self.text_widget.insert('end', f"  pNN50: {stats['pnn50']:.1f}%\n")
        
        self.text_widget.config(state='disabled')
    
    def export_to_csv(self, stats):
        if not stats:
            messagebox.showwarning("Предупреждение", "Нет данных для экспорта")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            try:
                df = pd.DataFrame([stats])
                df.to_csv(filename, index=False)
                messagebox.showinfo("Успех", f"Данные сохранены в {filename}")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при сохранении: {str(e)}")

class AnnotationWindow(tk.Toplevel):
    def __init__(self, parent, annotations=None):
        super().__init__(parent)
        self.title("Аннотации")
        self.geometry("600x400")
        self.transient(parent)
        self.grab_set()
        
        self.annotations = annotations or []
        
        # Список аннотаций
        list_frame = ttk.LabelFrame(self, text="Список аннотаций")
        list_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        self.annotation_list = tk.Listbox(list_frame, font=('Arial', 10))
        scrollbar = ttk.Scrollbar(list_frame, orient='vertical', command=self.annotation_list.yview)
        self.annotation_list.configure(yscrollcommand=scrollbar.set)
        
        self.annotation_list.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        # Кнопки управления
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', padx=10, pady=5)
        ttk.Button(btn_frame, text="Добавить", command=self.add_annotation).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Удалить", command=self.remove_annotation).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Редактировать", command=self.edit_annotation).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Закрыть", command=self.destroy).pack(side='right', padx=5)
        
        self.update_list()
    
    def update_list(self):
        self.annotation_list.delete(0, tk.END)
        for i, ann in enumerate(self.annotations):
            self.annotation_list.insert(tk.END, f"{i+1}. {ann['time']} - {ann['text']}")
    
    def add_annotation(self):
        text = simpledialog.askstring("Добавить аннотацию", "Введите текст аннотации:")
        if text:
            time_str = datetime.now().strftime("%H:%M:%S")
            self.annotations.append({'time': time_str, 'text': text})
            self.update_list()
    
    def remove_annotation(self):
        selection = self.annotation_list.curselection()
        if selection:
            idx = selection[0]
            del self.annotations[idx]
            self.update_list()
    
    def edit_annotation(self):
        selection = self.annotation_list.curselection()
        if selection:
            idx = selection[0]
            old_text = self.annotations[idx]['text']
            new_text = simpledialog.askstring("Редактировать аннотацию", 
                                            "Введите новый текст:", initialvalue=old_text)
            if new_text:
                self.annotations[idx]['text'] = new_text
                self.update_list()

class ExportWindow(tk.Toplevel):
    def __init__(self, parent, data_dict):
        super().__init__(parent)
        self.title("Экспорт данных")
        self.geometry("400x300")
        self.transient(parent)
        self.grab_set()
        
        self.data = data_dict
        
        # Выбор данных для экспорта
        ttk.Label(self, text="Выберите данные для экспорта:", font=('Arial', 12, 'bold')).pack(pady=10)
        
        self.export_avg = tk.BooleanVar(value=True)
        self.export_spectrum = tk.BooleanVar(value=True)
        self.export_rr_stats = tk.BooleanVar(value=True)
        self.export_annotations = tk.BooleanVar(value=True)
        
        ttk.Checkbutton(self, text="Усреднённый комплекс", variable=self.export_avg).pack(anchor='w', padx=20, pady=5)
        ttk.Checkbutton(self, text="Спектр", variable=self.export_spectrum).pack(anchor='w', padx=20, pady=5)
        ttk.Checkbutton(self, text="Статистика RR интервалов", variable=self.export_rr_stats).pack(anchor='w', padx=20, pady=5)
        ttk.Checkbutton(self, text="Аннотации", variable=self.export_annotations).pack(anchor='w', padx=20, pady=5)
        
        # Формат экспорта
        format_frame = ttk.LabelFrame(self, text="Формат экспорта")
        format_frame.pack(fill='x', padx=10, pady=10)
        
        self.export_format = tk.StringVar(value="csv")
        ttk.Radiobutton(format_frame, text="CSV", variable=self.export_format, value="csv").pack(anchor='w', padx=10, pady=2)
        ttk.Radiobutton(format_frame, text="JSON", variable=self.export_format, value="json").pack(anchor='w', padx=10, pady=2)
        ttk.Radiobutton(format_frame, text="MATLAB (.mat)", variable=self.export_format, value="mat").pack(anchor='w', padx=10, pady=2)
        
        # Кнопки
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill='x', padx=10, pady=10)
        ttk.Button(btn_frame, text="Экспорт", command=self.export_data).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Отмена", command=self.destroy).pack(side='right', padx=5)
    
    def export_data(self):
        filename = filedialog.asksaveasfilename(
            defaultextension=f".{self.export_format.get()}",
            filetypes=[
                ("CSV files", "*.csv"),
                ("JSON files", "*.json"),
                ("MATLAB files", "*.mat"),
                ("All files", "*.*")
            ]
        )
        
        if not filename:
            return
        
        try:
            export_data = {}
            
            if self.export_avg.get() and 'avg_complex' in self.data:
                export_data['avg_complex'] = {
                    'time': self.data['avg_complex']['time'].tolist(),
                    'amplitude': self.data['avg_complex']['amplitude'].tolist()
                }
            
            if self.export_spectrum.get() and 'spectrum' in self.data:
                export_data['spectrum'] = {
                    'frequency': self.data['spectrum']['frequency'].tolist(),
                    'amplitude': self.data['spectrum']['amplitude'].tolist()
                }
            
            if self.export_rr_stats.get() and 'rr_stats' in self.data:
                export_data['rr_statistics'] = self.data['rr_stats']
            
            if self.export_annotations.get() and 'annotations' in self.data:
                export_data['annotations'] = self.data['annotations']
            
            # Добавляем метаданные
            export_data['metadata'] = {
                'export_time': datetime.now().isoformat(),
                'sampling_rate': FS,
                'total_beats': self.data.get('total_beats', 0)
            }
            
            # Сохраняем в выбранном формате
            if self.export_format.get() == "csv":
                self.save_as_csv(filename, export_data)
            elif self.export_format.get() == "json":
                with open(filename, 'w', encoding='utf-8') as f:
                    json.dump(export_data, f, indent=2, ensure_ascii=False)
            elif self.export_format.get() == "mat":
                self.save_as_mat(filename, export_data)
            
            messagebox.showinfo("Успех", f"Данные экспортированы в {filename}")
            self.destroy()
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при экспорте: {str(e)}")
    
    def save_as_csv(self, filename, data):
        # Создаём отдельные CSV файлы для разных типов данных
        base_name = filename.replace('.csv', '')
        
        if 'avg_complex' in data:
            df_avg = pd.DataFrame({
                'time_ms': data['avg_complex']['time'],
                'amplitude_mv': data['avg_complex']['amplitude']
            })
            df_avg.to_csv(f"{base_name}_avg_complex.csv", index=False)
        
        if 'spectrum' in data:
            df_sp = pd.DataFrame({
                'frequency_hz': data['spectrum']['frequency'],
                'amplitude': data['spectrum']['amplitude']
            })
            df_sp.to_csv(f"{base_name}_spectrum.csv", index=False)
        
        if 'rr_statistics' in data:
            df_rr = pd.DataFrame([data['rr_statistics']])
            df_rr.to_csv(f"{base_name}_rr_statistics.csv", index=False)
        
        if 'annotations' in data:
            df_ann = pd.DataFrame(data['annotations'])
            df_ann.to_csv(f"{base_name}_annotations.csv", index=False)
    
    def save_as_mat(self, filename, data):
        try:
            from scipy.io import savemat
            mat_data = {}
            
            if 'avg_complex' in data:
                mat_data['avg_time'] = np.array(data['avg_complex']['time'])
                mat_data['avg_amplitude'] = np.array(data['avg_complex']['amplitude'])
            
            if 'spectrum' in data:
                mat_data['spectrum_freq'] = np.array(data['spectrum']['frequency'])
                mat_data['spectrum_amp'] = np.array(data['spectrum']['amplitude'])
            
            if 'rr_statistics' in data:
                for key, value in data['rr_statistics'].items():
                    mat_data[f'rr_{key}'] = value
            
            if 'metadata' in data:
                mat_data['metadata'] = data['metadata']
            
            savemat(filename, mat_data)
            
        except ImportError:
            messagebox.showerror("Ошибка", "Для экспорта в MATLAB формат требуется scipy.io")
            raise

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
        filemenu.add_command(label="Экспорт данных...", command=self.show_export_window)
        filemenu.add_separator()
        filemenu.add_command(label="Выход", command=lambda: on_close(self))
        menubar.add_cascade(label="Файл", menu=filemenu)
        
        # Меню Анализ
        analysismenu = tk.Menu(menubar, tearoff=0)
        analysismenu.add_command(label="Статистика RR интервалов", command=self.show_rr_statistics)
        analysismenu.add_command(label="Настройки анализа...", command=self.show_settings)
        analysismenu.add_separator()
        analysismenu.add_command(label="Сохранить настройки", command=self.save_settings)
        menubar.add_cascade(label="Анализ", menu=analysismenu)
        
        # Меню Аннотации
        annotationmenu = tk.Menu(menubar, tearoff=0)
        annotationmenu.add_command(label="Управление аннотациями...", command=self.show_annotations)
        menubar.add_cascade(label="Аннотации", menu=annotationmenu)
        
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
        
        # --- Новые переменные для дополнительных функций ---
        self.bad_beats = set()       # множество индексов плохих ударов
        self.manual_beats = set()    # множество индексов ручных ударов
        self.annotations = []        # список аннотаций
        self.rr_stats = {}           # статистика RR интервалов
        self.export_data = {}        # данные для экспорта
        
        # --- Параметры фильтров ---
        self.bandpass_enabled = False
        self.lowpass_enabled = False
        self.highpass_enabled = False
        self.bandpass_low = 0.5
        self.bandpass_high = 100.0
        self.lowpass_cutoff = 100.0
        self.highpass_cutoff = 0.5

        # --- Разметка grid ---
        self.rowconfigure(0, weight=3)
        self.rowconfigure(1, weight=2)
        self.rowconfigure(2, weight=2)
        self.rowconfigure(3, weight=0)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        # --- Raw ECG и HR ---
        self.fig_raw, self.ax_raw = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_raw = FigureCanvasTkAgg(self.fig_raw, master=self)
        self.canvas_raw.get_tk_widget().grid(row=0, column=0, columnspan=2, sticky='nsew')
        self.hr_label = ttk.Label(self.canvas_raw.get_tk_widget(), text="ЧСС: -- bpm", background='#fff')
        self.hr_label.place(relx=0.01, rely=0.01)

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

        # --- Выбранный комплекс ---
        self.fig_sel, self.ax_sel = plt.subplots(figsize=(8,2))
        plt.tight_layout()
        self.canvas_sel = FigureCanvasTkAgg(self.fig_sel, master=self)
        self.canvas_sel.get_tk_widget().grid(row=2, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

        # Загружаем сохранённые настройки
        self.load_settings()

        # --- Панель кнопок ---
        ctrl = ttk.Frame(self)
        ctrl.grid(row=3, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        
        # Первая строка кнопок
        ctrl_row1 = ttk.Frame(ctrl)
        ctrl_row1.pack(fill='x', pady=2)
        ttk.Button(ctrl_row1, text='Notch 50 Hz', command=self.toggle_notch).pack(side='left', padx=5)
        ttk.Button(ctrl_row1, text='Полосовой фильтр', command=self.toggle_bandpass).pack(side='left', padx=5)
        ttk.Button(ctrl_row1, text='Низкочастотный', command=self.toggle_lowpass).pack(side='left', padx=5)
        ttk.Button(ctrl_row1, text='Высокочастотный', command=self.toggle_highpass).pack(side='left', padx=5)
        
        # Вторая строка кнопок
        ctrl_row2 = ttk.Frame(ctrl)
        ctrl_row2.pack(fill='x', pady=2)
        ttk.Button(ctrl_row2, text='Выделить регион', command=self.activate_region_selector).pack(side='left', padx=5)
        ttk.Button(ctrl_row2, text='Отметить комплекс', command=self.toggle_peak_selector).pack(side='left', padx=5)
        ttk.Button(ctrl_row2, text='Выбрать комплекс', command=self.toggle_select_complex).pack(side='left', padx=5)
        ttk.Button(ctrl_row2, text='Статистика RR', command=self.show_rr_statistics).pack(side='left', padx=5)
        
        # Третья строка кнопок
        ctrl_row3 = ttk.Frame(ctrl)
        ctrl_row3.pack(fill='x', pady=2)
        ttk.Button(ctrl_row3, text='Обновить', command=self.compute).pack(side='left', padx=5)
        ttk.Button(ctrl_row3, text='Экспорт', command=self.show_export_window).pack(side='left', padx=5)
        ttk.Button(ctrl_row3, text='Аннотации', command=self.show_annotations).pack(side='left', padx=5)
        ttk.Button(ctrl_row3, text='Copy Screenshot', command=self.copy_screenshot).pack(side='left', padx=5)

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

    def toggle_bandpass(self):
        self.bandpass_enabled = not self.bandpass_enabled
        if self.bandpass_enabled:
            self.lowpass_enabled = False
            self.highpass_enabled = False
        self.compute()

    def toggle_lowpass(self):
        self.lowpass_enabled = not self.lowpass_enabled
        if self.lowpass_enabled:
            self.bandpass_enabled = False
            self.highpass_enabled = False
        self.compute()

    def toggle_highpass(self):
        self.highpass_enabled = not self.highpass_enabled
        if self.highpass_enabled:
            self.bandpass_enabled = False
            self.lowpass_enabled = False
        self.compute()

    def show_rr_statistics(self):
        if self.rr_stats:
            StatisticsWindow(self, self.rr_stats)
        else:
            messagebox.showinfo("Информация", "Сначала загрузите данные и выполните анализ")

    def show_settings(self):
        SettingsWindow(self)

    def show_annotations(self):
        AnnotationWindow(self, self.annotations)

    def show_export_window(self):
        if self.export_data:
            ExportWindow(self, self.export_data)
        else:
            messagebox.showinfo("Информация", "Сначала загрузите данные и выполните анализ")

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
        
        # Применяем фильтры
        if self.notch_enabled:
            ecg = notch_filter(ecg, FS)
        if self.bandpass_enabled:
            ecg = butter_bandpass_filter(ecg, FS, self.bandpass_low, self.bandpass_high)
        if self.lowpass_enabled:
            ecg = butter_lowpass_filter(ecg, FS, self.lowpass_cutoff)
        if self.highpass_enabled:
            ecg = butter_highpass_filter(ecg, FS, self.highpass_cutoff)
        
        self.ecg = ecg

        self.r_peaks = detect_r_peaks(ecg, FS)
        # HR по годным пикам
        good = [r for i,r in enumerate(self.r_peaks)
                if i not in self.bad_beats and not any(a<=r/FS<=b for a,b in self.bad_regions)]
        if len(good)>1:
            rr = np.diff(np.array(good)/FS)
            hr = (60/rr).mean()
            self.hr_label.config(text=f"ЧСС: {hr:.0f} bpm")
        else:
            self.hr_label.config(text="ЧСС: -- bpm")

        # Вычисляем статистику RR интервалов
        if len(good) > 1:
            self.rr_stats = calculate_rr_statistics(np.array(good), FS)
        else:
            self.rr_stats = {}

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

        # Сохраняем данные для экспорта
        self.export_data = {
            'avg_complex': {'time': t_avg, 'amplitude': avg} if avg is not None else {},
            'spectrum': {'frequency': xf[mask], 'amplitude': amp[mask]} if avg is not None else {},
            'rr_stats': self.rr_stats,
            'annotations': self.annotations,
            'total_beats': len(self.r_peaks)
        }

        self.draw_selected()

    def draw_raw(self):
        self.ax_raw.cla()
        t = np.arange(len(self.ecg)) / FS
        self.ax_raw.plot(t, self.ecg, color='black', linewidth=0.8)
        self.ax_raw.plot(self.r_peaks/FS, self.ecg[self.r_peaks], 'ro', markersize=3)
        for patch in self.region_patches:
            self.ax_raw.add_patch(patch)
        
        # Показываем статус фильтров в заголовке
        filter_status = []
        if self.notch_enabled:
            filter_status.append("Notch")
        if self.bandpass_enabled:
            filter_status.append("Bandpass")
        if self.lowpass_enabled:
            filter_status.append("Lowpass")
        if self.highpass_enabled:
            filter_status.append("Highpass")
        
        if filter_status:
            self.ax_raw.set_title(f'ЭКГ сигнал (Фильтры: {", ".join(filter_status)})')
        else:
            self.ax_raw.set_title('ЭКГ сигнал (Без фильтров)')
        
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

    def save_settings(self):
        """Сохранение настроек в файл"""
        settings = {
            'notch_enabled': self.notch_enabled,
            'bandpass_enabled': self.bandpass_enabled,
            'lowpass_enabled': self.lowpass_enabled,
            'highpass_enabled': self.highpass_enabled,
            'bandpass_low': self.bandpass_low,
            'bandpass_high': self.bandpass_high,
            'lowpass_cutoff': self.lowpass_cutoff,
            'highpass_cutoff': self.highpass_cutoff,
            'annotations': self.annotations
        }
        try:
            with open('ecg_settings.json', 'w', encoding='utf-8') as f:
                json.dump(settings, f, indent=2, ensure_ascii=False)
        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось сохранить настройки: {str(e)}")

    def load_settings(self):
        """Загрузка настроек из файла"""
        try:
            if os.path.exists('ecg_settings.json'):
                with open('ecg_settings.json', 'r', encoding='utf-8') as f:
                    settings = json.load(f)
                
                self.notch_enabled = settings.get('notch_enabled', True)
                self.bandpass_enabled = settings.get('bandpass_enabled', False)
                self.lowpass_enabled = settings.get('lowpass_enabled', False)
                self.highpass_enabled = settings.get('highpass_enabled', False)
                self.bandpass_low = settings.get('bandpass_low', 0.5)
                self.bandpass_high = settings.get('bandpass_high', 100.0)
                self.lowpass_cutoff = settings.get('lowpass_cutoff', 100.0)
                self.highpass_cutoff = settings.get('highpass_cutoff', 0.5)
                self.annotations = settings.get('annotations', [])
        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось загрузить настройки: {str(e)}")

if __name__ == '__main__':
    app = ECGApp()
    app.mainloop()




