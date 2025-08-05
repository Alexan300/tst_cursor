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
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Открыть запись...", command=self.open_file)
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

        # --- Панель кнопок ---
        ctrl = ttk.Frame(self)
        ctrl.grid(row=3, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        ttk.Button(ctrl, text='Notch 50 Hz', command=self.toggle_notch).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Выделить регион', command=self.activate_region_selector).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Отметить комплекс', command=self.toggle_peak_selector).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Выбрать комплекс', command=self.toggle_select_complex).pack(side='left', padx=5)
        ttk.Button(ctrl, text='Обновить', command=self.compute).pack(side='left', padx=5)
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




