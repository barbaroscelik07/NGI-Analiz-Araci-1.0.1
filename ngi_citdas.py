"""
NGI Cascade Impactor Analysis Tool — PyQt6
CITDAS-validated · Ph.Eur 2.9.18 / USP <601>
"""
import sys, os, math, datetime
import numpy as np
from scipy.stats import norm
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QSplitter,
    QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QLineEdit, QPushButton, QCheckBox, QComboBox,
    QRadioButton, QButtonGroup, QTabWidget, QScrollArea,
    QFrame, QFileDialog, QMessageBox, QDialog,
    QListWidget, QListWidgetItem, QGroupBox, QSizePolicy,
    QSpacerItem, QProgressBar, QDialogButtonBox
)
from PyQt6.QtCore import (Qt, QThread, pyqtSignal, QTimer, QSize)
from PyQt6.QtGui import (QColor, QPalette, QFont, QIcon, QPixmap)
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

# ─── resource_path ────────────────────────────────────────────────────────────
def resource_path(rel):
    base = getattr(sys, '_MEIPASS',
        os.path.dirname(os.path.abspath(
            sys.executable if getattr(sys,'frozen',False) else __file__)))
    return os.path.join(base, rel)

# ═══════════════════════════════════════════════════════════════════════════════
# HESAPLAMA MOTORU (değişmedi)
# ═══════════════════════════════════════════════════════════════════════════════
NGI_CUTOFFS = {
    15: {"Device":999,"Throat":999,"Presep":999,"S1":14.10,
         "S2":8.61,"S3":5.39,"S4":3.30,"S5":2.08,"S6":1.36,"S7":0.98,"MOC":0.54},
    30: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":11.719,"S3":6.395,"S4":3.988,"S5":2.299,"S6":1.357,"S7":0.834,"MOC":0.541},
    40: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":10.033,"S3":5.507,"S4":3.454,"S5":2.008,"S6":1.165,"S7":0.701,"MOC":0.446},
    60: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":8.06,"S3":4.46,"S4":2.82,"S5":1.66,"S6":0.94,"S7":0.55,"MOC":0.34},
    75: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":7.145,"S3":3.971,"S4":2.522,"S5":1.495,"S6":0.835,"S7":0.481,"MOC":0.293},
}
ALL_KEYS      = ["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]
ISM_STAGES    = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
GRAPH_STAGES  = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
DISP_STAGES   = ["Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]
RUNS_PER_SER  = 3
EXCL_FLOWS    = {15}
CP = ["#2E75B6","#ED7D31","#70AD47","#E84040","#7030A0",
      "#00B0F0","#D4A000","#C00000","#00B050","#FF69B4",
      "#FF8C00","#4169E1","#DC143C","#228B22","#8B008B"]

def calc_run(masses, flow, lo=15, hi=85, delivered_tp=False):
    co   = NGI_CUTOFFS[flow]
    excl = flow in EXCL_FLOWS
    ism  = sum(masses.get(s,0) for s in ISM_STAGES)
    metered = masses.get("Throat",0) + masses.get("Presep",0) + ism
    if ism <= 0:
        return {"error":"no_data","metered":metered,"delivered":ism,"masses":masses}
    cum = []
    for s in ALL_KEYS:
        u = 0.0
        if s in ISM_STAGES and co.get(s,999) < 900:
            if excl:
                u = sum(masses.get(x,0) for x in ISM_STAGES
                        if co.get(x,999) < co.get(s,999)) / ism * 100
            else:
                u = sum(masses.get(x,0) for x in ISM_STAGES
                        if co.get(x,999) <= co.get(s,999)) / ism * 100
        cum.append({"stage":s,"d50":co.get(s,999),"mass":masses.get(s,0),"u_pct":u})
    valid = [r for r in cum if r["stage"] in ISM_STAGES
             and co.get(r["stage"],999) < 900 and lo < r["u_pct"] < hi]
    delivered = (masses.get("Throat",0)+masses.get("Presep",0)+ism) if delivered_tp else ism
    res = {"metered":metered,"delivered":delivered,"cum_data":cum,
           "valid":valid,"masses":masses,"flow":flow}
    if len(valid) < 2:
        res["error"] = "insufficient"; res["n"] = len(valid); return res
    pts_all = sorted([(co[r["stage"]], r["u_pct"]) for r in cum
                      if r["stage"] in ISM_STAGES and co.get(r["stage"],999)<900],
                     key=lambda p: p[0])
    x = np.array([math.log10(v["d50"]) for v in valid])
    y = np.array([norm.ppf(v["u_pct"]/100) for v in valid])
    b = np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
    a = y.mean()-b*x.mean()
    yp = a+b*x; ss_r=np.sum((y-yp)**2); ss_t=np.sum((y-y.mean())**2)
    r2 = 1-ss_r/ss_t if (len(valid)>2 and ss_t>0) else 1.0
    mmad = 10**(-a/b)
    for i in range(len(pts_all)-1):
        d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
        if u1<=50<=u2 and u2>u1:
            t=(50-u1)/(u2-u1)
            mmad=10**(math.log10(d1)+t*(math.log10(d2)-math.log10(d1))); break
    def get_d(tu,tz):
        for i in range(len(pts_all)-1):
            d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
            if u1<=tu<=u2 and u2>u1:
                if any(lo<u<hi for _,u in [pts_all[i],pts_all[i+1]]):
                    t=(tu-u1)/(u2-u1)
                    return 10**(math.log10(d1)+t*(math.log10(d2)-math.log10(d1)))
        return 10**((tz-a)/b)
    d84=get_d(84.13,1.0); d16=get_d(15.87,-1.0)
    gsd=math.sqrt(d84/d16) if (d84 and d16 and d16>0) else 10**(1/b)
    def _interp_at(d_cut):
        """Cumulative % at given cutoff diameter via log interpolation"""
        for i in range(len(pts_all)-1):
            d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
            if d1<=d_cut<=d2 and d2>d1:
                t=(math.log10(d_cut)-math.log10(d1))/(math.log10(d2)-math.log10(d1))
                return u1+t*(u2-u1)
        return norm.cdf(a+b*math.log10(d_cut))*100

    d5u  = _interp_at(5.0)
    d3u  = _interp_at(3.0)
    d15u = _interp_at(1.5)
    fpd   = d5u /100*ism; fpf   = fpd /metered*100 if metered>0 else 0
    fpd3  = d3u /100*ism; fpf3  = fpd3/metered*100 if metered>0 else 0
    fpd15 = d15u/100*ism; fpf15 = fpd15/metered*100 if metered>0 else 0
    # D10 / D50 / D90 — kümülatif eğriden log-interp
    def _dv_at(target_pct):
        """Diameter at target cumulative % (log interpolation on pts_all)"""
        for i in range(len(pts_all)-1):
            d1,u1 = pts_all[i]; d2,u2 = pts_all[i+1]
            if u1 <= target_pct <= u2 and u2 > u1:
                t = (target_pct - u1) / (u2 - u1)
                return 10**(math.log10(d1) + t*(math.log10(d2)-math.log10(d1)))
        # Fallback: regresyon
        try:
            return 10**((norm.ppf(target_pct/100) - a) / b)
        except: return None
    d10 = _dv_at(10.0)
    d50 = _dv_at(50.0)  # MMAD ile aynı
    d90 = _dv_at(90.0)
    res.update({"n":len(valid),"a":a,"b":b,"slope":b,"intercept":a+5,"r2":r2,
                "mmad":mmad,"gsd":gsd,
                "fpd":fpd,"fpf":fpf,
                "fpd3":fpd3,"fpf3":fpf3,
                "fpd15":fpd15,"fpf15":fpf15,
                "d10":d10,"d50":d50,"d90":d90,
                "x_reg":x,"y_reg":y})
    return res

def calc_series_avg(runs):
    valid=[r for r in runs if "error" not in r]
    if not valid: return None
    avg_masses={}
    for s in ALL_KEYS:
        vals=[r["masses"].get(s,0) for r in valid]
        avg_masses[s]=float(np.mean(vals))
    params={}
    for p in ["mmad","gsd","fpd","fpf","fpd3","fpf3","fpd15","fpf15","d10","d50","d90","metered","delivered","slope","intercept","r2"]:
        vals=[r[p] for r in valid if p in r]
        if vals:
            m=float(np.mean(vals)); sd=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            params[p]=(m,sd,sd/m*100 if m else 0.0)
    return {"avg_masses":avg_masses,"params":params,"n_valid":len(valid)}

def calc_f2(ref_m, test_m, co):
    stages=[s for s in GRAPH_STAGES if co.get(s,999)<900]
    diffs=[]
    for s in stages:
        r=ref_m.get(s,0); t=test_m.get(s,0)
        if r>0: diffs.append(((t-r)/r*100)**2)
    if not diffs: return None
    return 50*math.log10(100/math.sqrt(1+np.mean(diffs)))

def fmt_num(v, decimals=4, dec_sep=","):
    if v is None: return "-"
    if isinstance(v,int): return str(v)
    try: return f"{v:.{decimals}f}".replace(".",dec_sep)
    except: return str(v)

def parse_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        raw = f.read()
    lines = [l.strip() for l in raw.splitlines() if l.strip()]
    if not lines: return None, None, ["Dosya bos"]
    sep = ";" if raw.count(";") > raw.count(",") else ","
    hdr_idx = None
    for i, line in enumerate(lines):
        cols = [c.strip().lower() for c in line.split(sep)]
        if "seri" in cols and ("run" in cols or "flow" in cols):
            hdr_idx = i; break
    if hdr_idx is None:
        return None, None, ["Header satiri bulunamadi. 'Seri', 'Run', 'Flow' kolonlari olmali."]
    hdr = [h.strip().lower().replace(" ","_") for h in lines[hdr_idx].split(sep)]
    lines = lines[hdr_idx+1:]
    missing = [c for c in ["seri","run","flow","throat","s1"] if c not in hdr]
    if missing: return None, None, [f"Eksik kolon: {missing}"]
    def to_float(s):
        if s is None or str(s).strip()=="": return 0.0
        v=str(s).strip()
        if v.count(",")==1 and "." not in v: v=v.replace(",",".")
        elif v.count(",")>1: v=v.replace(",","")
        try: return float(v)
        except: return 0.0
    ref_col="referans" in hdr
    series_dict={}; flow_vals=set(); warnings=[]
    for line in lines:
        if not line.replace(sep,"").strip(): continue
        parts=line.split(sep)
        row={hdr[i]:parts[i].strip() if i<len(parts) else "" for i in range(len(hdr))}
        seri=row.get("seri","").strip()
        if not seri: continue
        try:
            run_no=int(float(row.get("run","1") or "1"))
            flow_v=int(float(row.get("flow","60") or "60"))
        except: continue
        flow_vals.add(flow_v)
        is_ref=False
        if ref_col:
            rv=row.get("referans","0").strip()
            is_ref=rv in ("1","1.0","evet","yes","true","referans","ref")
        masses={"Device":0.0,"Throat":to_float(row.get("throat")),
                "Presep":to_float(row.get("presep")),
                "S1":to_float(row.get("s1")),"S2":to_float(row.get("s2")),
                "S3":to_float(row.get("s3")),"S4":to_float(row.get("s4")),
                "S5":to_float(row.get("s5")),"S6":to_float(row.get("s6")),
                "S7":to_float(row.get("s7")),"MOC":to_float(row.get("moc"))}
        if seri not in series_dict:
            series_dict[seri]={"runs":[],"ref":is_ref,"flow":flow_v}
        if is_ref: series_dict[seri]["ref"]=True
        if len(series_dict[seri]["runs"])<3:
            series_dict[seri]["runs"].append({"run_no":run_no,"masses":masses})
        else:
            key="csv_4runs__"+seri
            if key not in warnings: warnings.append(key)
    if len(flow_vals)>1: warnings.append("csv_err_flow")
    flow=sorted(flow_vals)[0] if flow_vals else 60
    return series_dict, flow, warnings

# ─── Çeviri ──────────────────────────────────────────────────────────────────
L = {
"TR":{
 "title":"NGI Impaktor Analiz Araci",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | CITDAS Dogrulamali",
 "lang_btn":"English","product":"Urun Adi","batch":"Lot No.",
 "operator":"Analist","date":"Tarih","flow_rate":"Akis Hizi",
 "add_series":"+ Seri Ekle","del_series":"Seri Sil",
 "calculate":"Hesapla","clear":"Temizle","export_pdf":"PDF Rapor",
 "load_csv":"CSV Yukle",
 "tab_scatter":"Nokta Grafik","tab_plot":"Log-Probit",
 "tab_dist":"Dagilim","tab_summary":"Ozet","tab_compare":"Karsilastirma",
 "series":"Seri","run":"Run","paste_btn":"Yapistir",
 "mean":"Ort.","sd":"SD","rsd":"RSD%","accept":"Kabul",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Yetersiz nokta","status_ready":"Hazir.",
 "status_done":"Hesaplama tamamlandi.","status_calc":"Hesaplaniyor...",
 "valid_range":"Gecerlilik (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"Bu seri REFERANS","ref_label":"REFERANS",
 "limit_label":"Limit Tipi","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manuel (%)","lim_pct":"Limit %",
 "f2_label":"f2 Benzerlik Faktoru","f2_pass":">=50 Benzer","f2_fail":"<50 Farkli",
 "outside_warn":"UYARI: Limit disi noktalar","no_ref":"Referans secilmedi",
 "ddu_label":"DDU Analizi","rsd_limit":"RSD Kabul (%)",
 "fp_dose":"FPD (mg)","fp_frac":"FPF (%)","slope_lbl":"Slope",
 "int_lbl":"Intercept","r2_lbl":"R2","n_lbl":"n",
 "param":"Parametre","stage":"Stage","mass_mg":"Kutle (mg)",
 "cum_mass":"Kum. Kutle","cum_pct":"Kum. %","valid_pt":"Gecerli",
 "probit_z":"Probit z","trend_mmad":"MMAD Trendi","trend_gsd":"GSD Trendi",
 "csv_ref_ask":"Referans seri secin","csv_ref_none":"Referans yok (atla)",
 "csv_loaded":"CSV yuklendi: {n} seri, {r} run",
 "csv_err_flow":"Uyari: Farkli flow hizlari tespit edildi",
 "csv_4runs":"Uyari: {s} serisi 3den fazla run iceriyor, ilk 3 alindi",
 "delivered_tp":"Delivered = T+P+ISM","lp_avg_only":"Sadece Seri Ortalamalari",
 "dec_sep":",",
},
"EN":{
 "title":"NGI Cascade Impactor Analysis",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | CITDAS Validated",
 "lang_btn":"Turkce","product":"Product","batch":"Batch No.",
 "operator":"Analyst","date":"Date","flow_rate":"Flow Rate",
 "add_series":"+ Add Series","del_series":"Del Series",
 "calculate":"Calculate","clear":"Clear","export_pdf":"PDF Report",
 "load_csv":"Load CSV",
 "tab_scatter":"Scatter Plot","tab_plot":"Log-Probit",
 "tab_dist":"Distribution","tab_summary":"Summary","tab_compare":"Compare",
 "series":"Series","run":"Run","paste_btn":"Paste",
 "mean":"Mean","sd":"SD","rsd":"RSD%","accept":"Accept",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Insufficient pts","status_ready":"Ready.",
 "status_done":"Calculation complete.","status_calc":"Calculating...",
 "valid_range":"Valid Range (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"This series is REFERENCE","ref_label":"REFERENCE",
 "limit_label":"Limit Type","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manual (%)","lim_pct":"Limit %",
 "f2_label":"f2 Similarity Factor","f2_pass":">=50 Similar","f2_fail":"<50 Different",
 "outside_warn":"WARNING: Points outside limits","no_ref":"No reference selected",
 "ddu_label":"DDU Analysis","rsd_limit":"RSD Accept (%)",
 "fp_dose":"FPD (mg)","fp_frac":"FPF (%)","slope_lbl":"Slope",
 "int_lbl":"Intercept","r2_lbl":"R2","n_lbl":"n",
 "param":"Parameter","stage":"Stage","mass_mg":"Mass (mg)",
 "cum_mass":"Cum. Mass","cum_pct":"Cum. %","valid_pt":"Valid",
 "probit_z":"Probit z","trend_mmad":"MMAD Trend","trend_gsd":"GSD Trend",
 "csv_ref_ask":"Select reference series","csv_ref_none":"No reference (skip)",
 "csv_loaded":"CSV loaded: {n} series, {r} runs",
 "csv_err_flow":"Warning: Multiple flow rates detected",
 "csv_4runs":"Warning: {s} has more than 3 runs, first 3 taken",
 "delivered_tp":"Delivered = T+P+ISM","lp_avg_only":"Series Averages Only",
 "dec_sep":".",
}}

# ─── Stil sabitleri ───────────────────────────────────────────────────────────
NAVY  = "#002D62"
NAVY2 = "#1a2e4a"
GOLD  = "#FFC600"
BG    = "#0e1219"
BG2   = "#141824"
BG3   = "#1c2336"
TXT   = "#e0eaf8"
TXT2  = "#7090b0"
GREEN = "#1a5a1a"
RED   = "#5a1a1a"

STYLE = f"""
QMainWindow, QDialog {{
    background: {BG};
}}
QWidget {{
    background: {BG};
    color: {TXT};
    font-family: 'Segoe UI', 'Arial';
    font-size: 13px;
}}
QLabel {{
    color: {TXT};
    background: transparent;
}}
QLineEdit {{
    background: {BG3};
    border: 1px solid #2a4060;
    border-radius: 4px;
    padding: 3px 6px;
    color: {TXT};
    selection-background-color: {NAVY2};
}}
QLineEdit:focus {{
    border: 1px solid {GOLD};
}}
QPushButton {{
    background: {NAVY2};
    border: 1px solid #2a4060;
    border-radius: 5px;
    padding: 5px 12px;
    color: {TXT};
    font-weight: bold;
}}
QPushButton:hover {{
    background: #253a5e;
    border-color: {GOLD};
}}
QPushButton:pressed {{
    background: #1a2a4a;
}}
QPushButton#btn_calc {{
    background: {GREEN};
    border-color: #2a8a2a;
    color: white;
}}
QPushButton#btn_calc:hover {{ background: #2a8a2a; }}
QPushButton#btn_pdf {{
    background: #3a1a5a;
    border-color: #6a20a0;
    color: white;
}}
QPushButton#btn_pdf:hover {{ background: #5a2a8a; }}
QPushButton#btn_csv {{
    background: #0a4a3a;
    border-color: #0a8a6a;
    color: white;
}}
QPushButton#btn_csv:hover {{ background: #0a7a5a; }}
QPushButton#btn_del {{
    background: {RED};
    border-color: #8a2020;
}}
QPushButton#btn_del:hover {{ background: #8a2020; }}
QPushButton#btn_clr {{
    background: #3a3a1a;
    border-color: #6a6a20;
}}
QComboBox {{
    background: {BG3};
    border: 1px solid #2a4060;
    border-radius: 4px;
    padding: 3px 6px;
    color: {TXT};
    min-width: 70px;
}}
QComboBox:focus {{ border-color: {GOLD}; }}
QComboBox QAbstractItemView {{
    background: {BG2};
    border: 1px solid #2a4060;
    color: {TXT};
    selection-background-color: {NAVY2};
}}
QComboBox::drop-down {{
    border: none;
    width: 20px;
}}
QCheckBox {{
    color: {TXT2};
    spacing: 6px;
}}
QCheckBox::indicator {{
    width: 14px; height: 14px;
    border: 1px solid #2a4060;
    border-radius: 3px;
    background: {BG3};
}}
QCheckBox::indicator:checked {{
    background: {NAVY2};
    border-color: {GOLD};
}}
QRadioButton {{
    color: {TXT2};
    spacing: 6px;
}}
QRadioButton::indicator {{
    width: 13px; height: 13px;
    border: 1px solid #2a4060;
    border-radius: 7px;
    background: {BG3};
}}
QRadioButton::indicator:checked {{
    background: {GOLD};
    border-color: {GOLD};
}}
QTabWidget::pane {{
    border: 1px solid #2a4060;
    background: {BG};
}}
QTabBar::tab {{
    background: {BG3};
    color: {TXT2};
    border: 1px solid #2a4060;
    border-bottom: none;
    padding: 8px 16px;
    margin-right: 2px;
    border-radius: 4px 4px 0 0;
    font-size: 13px;
}}
QTabBar::tab:selected {{
    background: #2E75B6;
    color: white;
    font-weight: bold;
}}
QTabBar::tab:hover:!selected {{
    background: #1a2e4a;
    color: {TXT};
}}
QScrollArea {{
    border: none;
    background: {BG};
}}
QScrollBar:vertical {{
    background: {BG2};
    width: 8px;
    border-radius: 4px;
}}
QScrollBar::handle:vertical {{
    background: #2a4060;
    border-radius: 4px;
    min-height: 20px;
}}
QScrollBar::handle:vertical:hover {{
    background: #3a6090;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0px;
}}
QScrollBar:horizontal {{
    background: {BG2};
    height: 8px;
    border-radius: 4px;
}}
QScrollBar::handle:horizontal {{
    background: #2a4060;
    border-radius: 4px;
    min-width: 20px;
}}
QScrollBar::handle:horizontal:hover {{
    background: #3a6090;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0px;
}}
QFrame[frameShape="4"], QFrame[frameShape="5"] {{
    color: #2a4060;
}}
QGroupBox {{
    border: 1px solid #2a4060;
    border-radius: 6px;
    margin-top: 8px;
    padding-top: 6px;
    color: {TXT2};
    font-size: 10px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 8px;
    padding: 0 4px;
}}
QSplitter::handle {{
    background: #2a4060;
    width: 2px;
}}
"""

# ═══════════════════════════════════════════════════════════════════════════════
# HESAPLAMA THREAD
# ═══════════════════════════════════════════════════════════════════════════════
class CalcThread(QThread):
    done = pyqtSignal(list)
    error = pyqtSignal(str)

    def __init__(self, series_inputs, flow, lo, hi, delivered_tp):
        super().__init__()
        self.series_inputs = series_inputs
        self.flow = flow; self.lo = lo; self.hi = hi
        self.delivered_tp = delivered_tp

    def run(self):
        try:
            results = []
            for si in self.series_inputs:
                runs = []
                for ri, masses in enumerate(si["masses_list"]):
                    r = calc_run(masses, self.flow, self.lo, self.hi, self.delivered_tp)
                    r["run_no"] = ri + 1
                    runs.append(r)
                avg = calc_series_avg(runs)
                results.append({
                    "name": si["name"], "color": si["color"],
                    "runs": runs, "avg": avg, "is_ref": si["is_ref"]
                })
            self.done.emit(results)
        except Exception as e:
            import traceback
            self.error.emit(traceback.format_exc())

# ═══════════════════════════════════════════════════════════════════════════════
# SERİ PANELİ
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# VERİ GİRİŞ POPUP
# ═══════════════════════════════════════════════════════════════════════════════
class DataEntryDialog(QDialog):
    def __init__(self, series_panel, T, color, parent=None):
        super().__init__(parent)
        self.sp = series_panel; self.T = T; self.color = color
        self.setWindowTitle(f"Veri Girisi - {series_panel.name_edit.text()}")
        self.setMinimumSize(600, 500)
        self.setStyleSheet(STYLE)
        self._build()
        self._load_existing()

    def _build(self):
        vl = QVBoxLayout(self)
        vl.setSpacing(6); vl.setContentsMargins(10,10,10,10)
        hdr = QLabel(f"  {self.sp.name_edit.text()} — Stage Kütleleri (mg)")
        hdr.setStyleSheet(f"color:{self.color};font-weight:bold;font-size:14px;"
            f"background:{BG3};border-left:4px solid {self.color};padding:6px;border-radius:4px;")
        vl.addWidget(hdr)
        pr = QHBoxLayout()
        lp = QLabel("Excel/Word'den kopyalayip yapistiriniz:")
        lp.setStyleSheet("color:#7090b0;font-size:11px;")
        pr.addWidget(lp); pr.addStretch()
        pb = QPushButton("Yapistir (Tum Sütunlar)")
        pb.setFixedSize(180,28); pb.setStyleSheet("font-size:11px;")
        pb.clicked.connect(self._paste_all); pr.addWidget(pb)
        vl.addLayout(pr)
        sc = QScrollArea(); sc.setWidgetResizable(True)
        gw = QWidget(); grid = QGridLayout(gw)
        grid.setSpacing(3); grid.setContentsMargins(4,4,4,4)
        for ci, txt in enumerate(["Stage","Run 1","Run 2","Run 3"]):
            l = QLabel(txt); l.setAlignment(Qt.AlignmentFlag.AlignCenter)
            l.setStyleSheet(f"color:{self.color if ci>0 else '#7090b0'};"
                f"font-weight:bold;font-size:12px;background:{BG3};border-radius:3px;padding:4px;")
            grid.addWidget(l, 0, ci)
        self.entries = [{} for _ in range(RUNS_PER_SER)]
        for si, s in enumerate(DISP_STAGES):
            lc = "#FFD700" if s=="Presep" else "#aac8e8"
            lbl = QLabel(s)
            lbl.setStyleSheet(f"color:{lc};font-size:13px;font-weight:bold;"
                "min-width:60px;background:transparent;")
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lbl, si+1, 0)
            for ri in range(RUNS_PER_SER):
                e = QLineEdit("0.0000")
                e.setAlignment(Qt.AlignmentFlag.AlignCenter)
                e.setFixedHeight(30)
                e.setStyleSheet("font-size:12px;padding:2px 6px;")
                e.focusInEvent = lambda ev, _e=e: (_e.selectAll() or
                    QLineEdit.focusInEvent(_e, ev))
                grid.addWidget(e, si+1, ri+1)
                self.entries[ri][s] = e
        sc.setWidget(gw); vl.addWidget(sc, 1)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel)
        btns.button(QDialogButtonBox.StandardButton.Ok).setText("Kaydet")
        btns.button(QDialogButtonBox.StandardButton.Cancel).setText("Iptal")
        btns.accepted.connect(self._save); btns.rejected.connect(self.reject)
        vl.addWidget(btns)

    def _load_existing(self):
        for ri in range(RUNS_PER_SER):
            for s in DISP_STAGES:
                self.entries[ri][s].setText(self.sp.entries[ri][s].text())

    def _save(self):
        for ri in range(RUNS_PER_SER):
            for s in DISP_STAGES:
                self.sp.entries[ri][s].setText(self.entries[ri][s].text())
        self.accept()

    def _paste_all(self):
        import re
        text = QApplication.clipboard().text()
        if not text.strip(): return
        lines = [l.strip() for l in text.strip().splitlines() if l.strip()]
        def split_line(ln):
            if '\t' in ln: return [t.strip() for t in ln.split('\t')]
            return [t.strip() for t in re.split(r'\s{2,}|\s+', ln.strip())]
        def is_hdr(ln):
            p=split_line(ln)
            try: float(p[0].replace(',','.').replace(' ','')); return False
            except: return True
        if lines and is_hdr(lines[0]): lines=lines[1:]
        all_s=["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]
        for ri, line in enumerate(lines[:RUNS_PER_SER]):
            tokens=split_line(line)
            try:
                vals=[float(t.replace(',','.').replace(' ','')) for t in tokens if t]
                if len(vals)>=10:
                    if len(vals)<11: vals=[0.0]+vals
                    for si,val in enumerate(vals[:11]):
                        s=all_s[si]
                        if s in self.entries[ri]:
                            self.entries[ri][s].setText(f"{val:.4f}")
            except: pass


# ═══════════════════════════════════════════════════════════════════════════════
# SERİ PANELİ — Kompakt satır
# ═══════════════════════════════════════════════════════════════════════════════
class SeriesPanel(QFrame):
    def __init__(self, idx, color, T, parent=None):
        super().__init__(parent)
        self.idx=idx; self.color=color; self.T=T
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setStyleSheet("""
            QFrame {
                background: rgba(255,255,255,0.04);
                border: 0.5px solid rgba(255,255,255,0.08);
                border-radius: 6px;
            }
            QFrame:hover {
                background: rgba(255,255,255,0.06);
                border: 0.5px solid rgba(255,198,0,0.2);
            }
        """)
        self.setFixedHeight(30)
        self.entries=[{s:QLineEdit("0.000") for s in DISP_STAGES}
                       for _ in range(RUNS_PER_SER)]
        self._build()

    def _build(self):
        hl = QHBoxLayout(self)
        hl.setSpacing(4); hl.setContentsMargins(6,3,6,3)
        # Renk çubuğu
        bar=QFrame(); bar.setFixedSize(3,22)
        bar.setStyleSheet(f"background:{self.color};border-radius:1px;")
        hl.addWidget(bar)
        # Görünürlük checkbox
        self.vis_check=QCheckBox(); self.vis_check.setChecked(True)
        self.vis_check.setToolTip("Grafikte göster/gizle")
        self.vis_check.setStyleSheet(f"""
            QCheckBox::indicator{{width:13px;height:13px;border:1.5px solid {self.color};
                border-radius:2px;background:transparent;}}
            QCheckBox::indicator:checked{{background:{self.color};}}
        """)
        hl.addWidget(self.vis_check)
        # İsim
        self.name_edit=QLineEdit(f"Seri {self.idx}")
        self.name_edit.setStyleSheet(f"font-weight:bold;font-size:12px;"
            f"background:transparent;border:none;color:#e0eaf8;min-width:90px;")
        self.name_edit.setFrame(False)
        hl.addWidget(self.name_edit,1)
        # Referans
        if self.idx==1:
            self.ref_check=QCheckBox("REF")
            self.ref_check.setStyleSheet(f"color:{GOLD};font-weight:bold;font-size:11px;")
            hl.addWidget(self.ref_check)
        else:
            self.ref_check=None
        # Düzenle
        self.edit_btn=QPushButton("✏")
        self.edit_btn.setFixedSize(22,22)
        self.edit_btn.setToolTip("Stage verilerini gir / Edit data")
        self.edit_btn.setStyleSheet("""
            QPushButton{background:transparent;border:none;
                color:#2a5a90;font-size:13px;border-radius:3px;}
            QPushButton:hover{color:#4a9ae8;background:rgba(20,50,100,0.3);}
        """)
        self.edit_btn.clicked.connect(self._open_edit)
        hl.addWidget(self.edit_btn)
        # Üç nokta menü (sil)
        self.del_btn=QPushButton("⋯")
        self.del_btn.setFixedSize(22,22)
        self.del_btn.setToolTip("Seriyi sil")
        self.del_btn.setStyleSheet("""
            QPushButton{background:transparent;border:none;
                color:#3a5070;font-size:15px;font-weight:bold;border-radius:3px;}
            QPushButton:hover{color:#FF6060;background:rgba(100,20,20,0.3);}
        """)
        hl.addWidget(self.del_btn)

    def _open_edit(self):
        dlg=DataEntryDialog(self, self.T, self.color, self)
        dlg.exec()

    def get_masses(self, ri):
        m={"Device":0.0}
        for s in DISP_STAGES:
            try: m[s]=float(self.entries[ri][s].text().replace(",","."))
            except: m[s]=0.0
        return m

    def set_masses(self, ri, masses):
        for s in DISP_STAGES:
            self.entries[ri][s].setText(fmt_num(masses.get(s,0.0),4,"."))

    def update_lang(self, T):
        self.T=T
        if self.ref_check:
            self.ref_check.setText("REF")


# ═══════════════════════════════════════════════════════════════════════════════
# HESAPLAMA THREAD
# ═══════════════════════════════════════════════════════════════════════════════
class CalcThread(QThread):
    done=pyqtSignal(list); error=pyqtSignal(str)
    def __init__(self,series_inputs,flow,lo,hi,delivered_tp):
        super().__init__()
        self.si=series_inputs; self.flow=flow
        self.lo=lo; self.hi=hi; self.dtp=delivered_tp
    def run(self):
        try:
            res=[]
            for si in self.si:
                runs=[]
                for ri,masses in enumerate(si["masses_list"]):
                    r=calc_run(masses,self.flow,self.lo,self.hi,self.dtp)
                    r["run_no"]=ri+1; runs.append(r)
                res.append({"name":si["name"],"color":si["color"],
                    "runs":runs,"avg":calc_series_avg(runs),"is_ref":si["is_ref"]})
            self.done.emit(res)
        except Exception:
            import traceback; self.error.emit(traceback.format_exc())


# ═══════════════════════════════════════════════════════════════════════════════
# ANA UYGULAMA — Seçenek C (Minimal & Kompakt)
# ═══════════════════════════════════════════════════════════════════════════════
class NGIApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.lang="TR"; self.T=L["TR"]
        self.all_series=[]; self.series_panels=[]
        self.flow=60; self.limit_type="ema"
        self.calc_thread=None
        self._sc_ser_checks = {}
        self._sec_series = None
        self._threads = []
        self._setup(); self._build_ui()
        self._add_series()

    def _setup(self):
        self.setWindowTitle(self.T["title"])
        self.resize(1540,980); self.setMinimumSize(1200,780)
        self.setStyleSheet(STYLE)
        ico=resource_path("icon.ico")
        if os.path.exists(ico):
            try: self.setWindowIcon(QIcon(ico))
            except: pass

    # ── UI ────────────────────────────────────────────────────────────────────
    def _build_ui(self):
        cw=QWidget(); self.setCentralWidget(cw)
        ml=QVBoxLayout(cw); ml.setSpacing(0); ml.setContentsMargins(0,0,0,0)
        ml.addWidget(self._build_header())
        spl=QSplitter(Qt.Orientation.Horizontal); spl.setHandleWidth(2)
        lw=QWidget(); lw.setFixedWidth(290)
        ll=QVBoxLayout(lw); ll.setContentsMargins(0,0,0,0); ll.setSpacing(0)
        self._build_left(ll); spl.addWidget(lw)
        rw=QWidget(); rl=QVBoxLayout(rw)
        rl.setContentsMargins(4,4,8,4); rl.setSpacing(0)
        self._build_right(rl); spl.addWidget(rw)
        spl.setSizes([290,1250]); ml.addWidget(spl,1)
        self.status_lbl=QLabel(self.T["status_ready"])
        self.status_lbl.setStyleSheet(
            f"color:{TXT2};font-size:11px;background:{BG};"
            f"padding:3px 10px;border-top:1px solid #1a2a40;")
        ml.addWidget(self.status_lbl)

    def _build_header(self):
        h=QWidget(); h.setFixedHeight(44)
        h.setStyleSheet(f"background:{NAVY};")
        hl=QHBoxLayout(h); hl.setContentsMargins(12,4,12,4)
        self.lbl_title=QLabel(self.T["title"])
        self.lbl_title.setStyleSheet(
            f"color:{GOLD};font-size:14px;font-weight:bold;background:transparent;")
        hl.addWidget(self.lbl_title)
        self.lbl_sub=QLabel(self.T["subtitle"])
        self.lbl_sub.setStyleSheet(
            "color:#6a9ec0;font-size:10px;background:transparent;")
        hl.addWidget(self.lbl_sub); hl.addStretch()
        self.btn_lang=QPushButton(self.T["lang_btn"])
        self.btn_lang.setFixedSize(76,26)
        self.btn_lang.setStyleSheet(
            "background:#001a40;border:1px solid #2a4060;"
            "color:#aac8e8;border-radius:3px;font-size:11px;")
        self.btn_lang.clicked.connect(self._toggle_lang)
        hl.addWidget(self.btn_lang)
        return h

    def _build_left(self, layout):
        layout.setSpacing(0); layout.setContentsMargins(0,0,0,0)
        # Tüm sol panel scroll
        scroll=QScrollArea(); scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet("background:#0d1117;border:none;")
        inner=QWidget(); inner.setStyleSheet("background:#0d1117;")
        vl=QVBoxLayout(inner); vl.setSpacing(6); vl.setContentsMargins(8,8,8,8)

        # ── Glassmorphic bölüm yardımcısı ────────────────────────────────────
        def glass_section(icon_char, title, expanded=True):
            outer=QFrame()
            outer.setStyleSheet("""
                QFrame {
                    background: rgba(255,255,255,0.04);
                    border: 0.5px solid rgba(255,255,255,0.09);
                    border-radius: 8px;
                }
            """)
            vbox=QVBoxLayout(outer); vbox.setSpacing(0); vbox.setContentsMargins(0,0,0,0)
            # Başlık
            hdr_btn=QPushButton()
            hdr_btn.setCheckable(True); hdr_btn.setChecked(expanded)
            hdr_btn.setStyleSheet(f"""
                QPushButton {{
                    background: rgba(255,255,255,0.0);
                    border: none;
                    border-bottom: 0.5px solid rgba(255,255,255,0.07);
                    border-radius: 8px 8px 0 0;
                    text-align: left;
                    padding: 8px 12px;
                    color: #c0d8f0;
                    font-size: 12px;
                    font-weight: bold;
                }}
                QPushButton:hover {{ background: rgba(255,255,255,0.04); }}
                QPushButton:checked {{
                    border-left: 2px solid {GOLD};
                    color: {GOLD};
                    background: rgba(0,45,98,0.25);
                }}
            """)
            hdr_inner=QHBoxLayout()
            hdr_inner.setContentsMargins(0,0,0,0); hdr_inner.setSpacing(8)
            ico=QLabel(icon_char)
            ico.setStyleSheet(f"color:{GOLD};font-size:13px;background:transparent;border:none;")
            hdr_inner.addWidget(ico)
            title_lbl=QLabel(title)
            title_lbl.setStyleSheet("color:#c0d8f0;font-size:12px;font-weight:bold;background:transparent;border:none;")
            hdr_inner.addWidget(title_lbl)
            hdr_inner.addStretch()
            arr=QLabel("▾")
            arr.setStyleSheet("color:#3a5070;font-size:10px;background:transparent;border:none;")
            hdr_inner.addWidget(arr)
            hdr_wgt=QWidget(); hdr_wgt.setLayout(hdr_inner)
            hdr_wgt.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
            hdr_btn_lay=QVBoxLayout(hdr_btn)
            hdr_btn_lay.setContentsMargins(0,0,0,0)
            hdr_btn_lay.addWidget(hdr_wgt)
            vbox.addWidget(hdr_btn)
            # Gövde
            body=QWidget()
            body.setStyleSheet("background:transparent;border:none;")
            body_lay=QVBoxLayout(body)
            body_lay.setSpacing(0); body_lay.setContentsMargins(0,0,0,0)
            body.setVisible(expanded)
            vbox.addWidget(body)
            def toggle(checked,bw=body,b=hdr_btn,a=arr,t=title,tl=title_lbl):
                bw.setVisible(checked)
                a.setText("▾" if checked else "▸")
            hdr_btn.toggled.connect(toggle)
            return {"outer":outer,"body":body_lay,"title_lbl":title_lbl}

        # ── Analiz Bilgileri ──────────────────────────────────────────────────
        sec_meta=glass_section("📋","Analiz Bilgileri")
        meta_w=QWidget(); meta_w.setStyleSheet("background:transparent;")
        ml2=QFormLayout(meta_w); ml2.setSpacing(5)
        ml2.setContentsMargins(10,6,10,10)
        ml2.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        for lbl_style in []:
            pass
        self.e_product=QLineEdit()
        self.e_batch=QLineEdit()
        self.e_operator=QLineEdit()
        self.e_date=QLineEdit(datetime.datetime.now().strftime("%d.%m.%Y"))
        for w in [self.e_product,self.e_batch,self.e_operator,self.e_date]:
            w.setStyleSheet("""
                QLineEdit{background:rgba(255,255,255,0.06);
                    border:0.5px solid rgba(255,255,255,0.12);
                    border-radius:4px;padding:4px 8px;color:#c0d8f0;}
                QLineEdit:focus{border:0.5px solid #FFC600;}
            """)
        ml2.addRow(self.T["product"], self.e_product)
        ml2.addRow(self.T["batch"],   self.e_batch)
        ml2.addRow(self.T["operator"],self.e_operator)
        ml2.addRow(self.T["date"],    self.e_date)
        # Form label stilleri
        for i in range(ml2.rowCount()):
            lbl_item=ml2.itemAt(i, QFormLayout.ItemRole.LabelRole)
            if lbl_item and lbl_item.widget():
                lbl_item.widget().setStyleSheet("color:#5a8ab0;font-size:11px;background:transparent;")
        sec_meta["body"].addWidget(meta_w)
        vl.addWidget(sec_meta["outer"])

        # ── Akış & Parametreler ───────────────────────────────────────────────
        sec_flow=glass_section("⚙","Akış & Parametreler")
        flow_w=QWidget(); flow_w.setStyleSheet("background:transparent;")
        fl=QVBoxLayout(flow_w); fl.setSpacing(5); fl.setContentsMargins(10,6,10,10)
        glass_input_style="""
            QLineEdit{background:rgba(255,255,255,0.06);
                border:0.5px solid rgba(255,255,255,0.12);
                border-radius:4px;padding:3px 6px;color:#c0d8f0;}
            QLineEdit:focus{border:0.5px solid #FFC600;}
        """
        glass_combo_style="""
            QComboBox{background:rgba(255,255,255,0.06);
                border:0.5px solid rgba(255,255,255,0.12);
                border-radius:4px;padding:3px 6px;color:#FFC600;font-weight:bold;}
            QComboBox:focus{border:0.5px solid #FFC600;}
            QComboBox QAbstractItemView{background:#1c2336;color:#c0d8f0;
                border:1px solid #2a4060;}
        """
        r1=QHBoxLayout(); r1.setSpacing(6)
        lf=QLabel(self.T["flow_rate"])
        lf.setStyleSheet(f"color:{GOLD};font-weight:bold;background:transparent;")
        r1.addWidget(lf)
        self.flow_combo=QComboBox()
        for fv in sorted(NGI_CUTOFFS.keys()): self.flow_combo.addItem(str(fv))
        self.flow_combo.setCurrentText("60"); self.flow_combo.setFixedWidth(72)
        self.flow_combo.setStyleSheet(glass_combo_style)
        self.flow_combo.currentTextChanged.connect(self._on_flow)
        r1.addWidget(self.flow_combo)
        lmin=QLabel("L/min"); lmin.setStyleSheet("color:#7090b0;background:transparent;")
        r1.addWidget(lmin); r1.addStretch()
        fl.addLayout(r1)
        r2=QHBoxLayout(); r2.setSpacing(4)
        lvr=QLabel(self.T["valid_range"]); lvr.setStyleSheet("color:#7090b0;background:transparent;")
        r2.addWidget(lvr)
        self.e_lo=QLineEdit("15"); self.e_lo.setFixedWidth(44)
        self.e_lo.setStyleSheet(glass_input_style)
        r2.addWidget(self.e_lo)
        dash=QLabel("–"); dash.setStyleSheet("color:#5a8ab0;background:transparent;")
        r2.addWidget(dash)
        self.e_hi=QLineEdit("85"); self.e_hi.setFixedWidth(44)
        self.e_hi.setStyleSheet(glass_input_style)
        r2.addWidget(self.e_hi); r2.addStretch()
        fl.addLayout(r2)
        r3=QHBoxLayout(); r3.setSpacing(6)
        lrsd=QLabel(self.T["rsd_limit"]); lrsd.setStyleSheet("color:#7090b0;background:transparent;")
        r3.addWidget(lrsd)
        self.e_rsd=QLineEdit("5"); self.e_rsd.setFixedWidth(38)
        self.e_rsd.setStyleSheet(glass_input_style)
        r3.addWidget(self.e_rsd)
        lpct=QLabel("%"); lpct.setStyleSheet("color:#5a8ab0;background:transparent;")
        r3.addWidget(lpct)
        self.chk_delivered_tp=QCheckBox(self.T["delivered_tp"])
        self.chk_delivered_tp.setStyleSheet("font-size:10px;color:#7090b0;")
        r3.addWidget(self.chk_delivered_tp); r3.addStretch()
        fl.addLayout(r3)
        co_w=QWidget()
        co_w.setStyleSheet("background:rgba(0,0,0,0.2);border-radius:4px;")
        self.cbox_layout=QHBoxLayout(co_w)
        self.cbox_layout.setContentsMargins(4,2,4,2); self.cbox_layout.setSpacing(2)
        self.cbox=co_w; fl.addWidget(co_w)
        sec_flow["body"].addWidget(flow_w)
        vl.addWidget(sec_flow["outer"])
        self._refresh_cutoffs()

        # ── İşlemler ─────────────────────────────────────────────────────────
        sec_ops=glass_section("▶","İşlemler")
        ops_w=QWidget(); ops_w.setStyleSheet("background:transparent;")
        ol=QVBoxLayout(ops_w); ol.setSpacing(4); ol.setContentsMargins(10,6,10,10)
        glass_btn = lambda text,obj=None,color=None: self._make_glass_btn(text,obj,color)

        def make_btn(text, obj_name=None, color=None, height=30):
            b=QPushButton(text); b.setFixedHeight(height)
            base=color or "#1a2e4a"
            b.setStyleSheet(f"""
                QPushButton{{background:{base};border:0.5px solid rgba(255,255,255,0.12);
                    border-radius:5px;color:#c0d8f0;font-size:12px;font-weight:bold;
                    padding:0 8px;}}
                QPushButton:hover{{background:rgba(255,255,255,0.08);
                    border:0.5px solid rgba(255,198,0,0.4);}}
                QPushButton:pressed{{background:rgba(0,0,0,0.2);}}
            """)
            if obj_name: b.setObjectName(obj_name)
            return b

        row1=QHBoxLayout(); row1.setSpacing(5)
        self.btn_calc=make_btn(self.T["calculate"],"btn_calc","rgba(20,80,20,0.7)")
        self.btn_calc.clicked.connect(self._calculate)
        row1.addWidget(self.btn_calc)
        self.btn_clr=make_btn(self.T["clear"],None,"rgba(60,50,10,0.6)")
        self.btn_clr.setFixedWidth(74); self.btn_clr.clicked.connect(self._clear)
        row1.addWidget(self.btn_clr)
        ol.addLayout(row1)

        row2=QHBoxLayout(); row2.setSpacing(5)
        self.btn_pdf=make_btn(self.T["export_pdf"],None,"rgba(60,20,90,0.7)")
        self.btn_pdf.clicked.connect(self._export_pdf)
        row2.addWidget(self.btn_pdf)
        self.btn_csv=make_btn(self.T["load_csv"],None,"rgba(10,60,50,0.7)")
        self.btn_csv.clicked.connect(self._load_csv)
        row2.addWidget(self.btn_csv)
        ol.addLayout(row2)
        sec_ops["body"].addWidget(ops_w)
        vl.addWidget(sec_ops["outer"])

        # ── Seriler ───────────────────────────────────────────────────────────
        self._sec_series=glass_section("◈","Seriler (0)")
        ser_w=QWidget(); ser_w.setStyleSheet("background:transparent;")
        sv=QVBoxLayout(ser_w); sv.setSpacing(4); sv.setContentsMargins(10,6,10,10)

        # Seri kontrol butonları
        sc_row=QHBoxLayout(); sc_row.setSpacing(4)
        self.btn_add=make_btn("+ "+self.T["add_series"],None,"rgba(20,50,100,0.6)",26)
        self.btn_add.clicked.connect(self._add_series)
        sc_row.addWidget(self.btn_add)
        self.btn_sel_all=make_btn("✓ Tümü",None,"rgba(20,50,20,0.5)",26)
        self.btn_sel_all.setFixedWidth(68)
        self.btn_sel_all.clicked.connect(self._select_all_series)
        sc_row.addWidget(self.btn_sel_all)
        self.btn_sel_none=make_btn("✕ Gizle",None,"rgba(60,20,20,0.5)",26)
        self.btn_sel_none.setFixedWidth(68)
        self.btn_sel_none.clicked.connect(self._deselect_all_series)
        sc_row.addWidget(self.btn_sel_none)
        sv.addLayout(sc_row)

        # Seri listesi
        self.series_container=QWidget()
        self.series_container.setStyleSheet("background:transparent;")
        self.series_layout=QVBoxLayout(self.series_container)
        self.series_layout.setSpacing(2); self.series_layout.setContentsMargins(0,0,0,0)
        self.series_layout.addStretch()
        sv.addWidget(self.series_container)
        self._sec_series["body"].addWidget(ser_w)
        vl.addWidget(self._sec_series["outer"])
        vl.addStretch()
        scroll.setWidget(inner); layout.addWidget(scroll)

    def _update_series_count(self):
        n=len(self.series_panels)
        lbl="Seriler" if self.lang=="TR" else "Series"
        if self._sec_series and "title_lbl" in self._sec_series:
            self._sec_series["title_lbl"].setText(f"{lbl} ({n})")


    def _build_right(self, layout):
        self.tabs=QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.currentChanged.connect(self._on_tab_change)
        layout.addWidget(self.tabs,1)

        # Nokta Grafik sekmesi
        self.tab_scatter=QWidget()
        sl=QVBoxLayout(self.tab_scatter); sl.setContentsMargins(4,4,4,4); sl.setSpacing(4)
        # Kontrol paneli
        sc_ctrl=QFrame(); sc_ctrl.setFixedHeight(68)
        sc_ctrl.setStyleSheet(f"background:{BG3};border-radius:4px;border:1px solid #2a4060;")
        sc_vl=QVBoxLayout(sc_ctrl); sc_vl.setContentsMargins(8,4,8,4); sc_vl.setSpacing(4)
        # Satır 1: parametre butonları
        row1=QHBoxLayout(); row1.setSpacing(4)
        row1.addWidget(QLabel("Parametre:"))
        self.scatter_param_grp=QButtonGroup(self)
        self._sc_param_btns={}
        for key,lbl in [("fpd","FPD"),("fpf","FPF"),("mmad","MMAD"),("gsd","GSD"),
                          ("d10","D10"),("d50","D50"),("d90","D90")]:
            btn=QPushButton(lbl); btn.setCheckable(True)
            btn.setFixedHeight(26)
            btn.setStyleSheet("""
                QPushButton{background:#1a2a4a;border:1px solid #2a4060;border-radius:4px;
                    color:#aac8e8;font-size:12px;font-weight:bold;padding:0 12px;}
                QPushButton:checked{background:#2E75B6;border-color:#4a95d6;color:white;}
                QPushButton:hover:!checked{background:#253a5e;}
            """)
            if key=="fpd": btn.setChecked(True)
            self.scatter_param_grp.addButton(btn)
            self._sc_param_btns[key]=btn
            row1.addWidget(btn)
        row1.addStretch()
        sc_vl.addLayout(row1)
        # Satır 2: alt ölçek butonları (FPD/FPF için)
        row2=QHBoxLayout(); row2.setSpacing(4)
        row2.addWidget(QLabel("Alt ölçek:"))
        self.scatter_scale_grp=QButtonGroup(self)
        self._sc_scale_btns={}
        for key,lbl in [("5","5 µm"),("3","3 µm"),("15","1,5 µm")]:
            btn=QPushButton(lbl); btn.setCheckable(True)
            btn.setFixedHeight(24)
            btn.setStyleSheet("""
                QPushButton{background:#141824;border:1px solid #2a4060;border-radius:4px;
                    color:#7090b0;font-size:11px;padding:0 10px;}
                QPushButton:checked{background:#1a3a6a;border-color:#2a5a9a;color:#aac8e8;}
                QPushButton:hover:!checked{background:#1a2a3a;}
            """)
            if key=="5": btn.setChecked(True)
            self.scatter_scale_grp.addButton(btn)
            self._sc_scale_btns[key]=btn
            row2.addWidget(btn)
        row2.addStretch()
        self.sc_row2_widget=QWidget()
        self.sc_row2_widget.setLayout(row2)
        sc_vl.addWidget(self.sc_row2_widget)
        sl.addWidget(sc_ctrl)
        # Grafik + seri listesi yan yana
        sc_body=QHBoxLayout(); sc_body.setSpacing(6)
        self.scatter_canvas=FigureCanvas(Figure(figsize=(8,5), facecolor=BG))
        sc_body.addWidget(self.scatter_canvas, 1)
        # Seri seçim paneli
        sc_ser_frame=QFrame()
        sc_ser_frame.setFixedWidth(180)
        sc_ser_frame.setStyleSheet(
            f"background:{BG3};border-radius:4px;border:1px solid #2a4060;")
        sc_ser_vl=QVBoxLayout(sc_ser_frame)
        sc_ser_vl.setContentsMargins(6,6,6,6); sc_ser_vl.setSpacing(4)
        sc_ser_title=QLabel("Seriler")
        sc_ser_title.setStyleSheet(
            f"color:{GOLD};font-weight:bold;font-size:12px;background:transparent;")
        sc_ser_vl.addWidget(sc_ser_title)
        sc_btn_row=QHBoxLayout(); sc_btn_row.setSpacing(4)
        sc_all=QPushButton("Tümü ✓"); sc_all.setFixedHeight(22)
        sc_all.setStyleSheet("font-size:10px;")
        sc_all.clicked.connect(self._scatter_sel_all)
        sc_none=QPushButton("Tümü ✕"); sc_none.setFixedHeight(22)
        sc_none.setStyleSheet("font-size:10px;background:#3a1a1a;")
        sc_none.clicked.connect(self._scatter_sel_none)
        sc_btn_row.addWidget(sc_all); sc_btn_row.addWidget(sc_none)
        sc_ser_vl.addLayout(sc_btn_row)
        # Seri listesi (scroll)
        sc_ser_scroll=QScrollArea(); sc_ser_scroll.setWidgetResizable(True)
        sc_ser_scroll.setStyleSheet("background:transparent;border:none;")
        self._sc_ser_container=QWidget()
        self._sc_ser_container.setStyleSheet("background:transparent;")
        self._sc_ser_layout=QVBoxLayout(self._sc_ser_container)
        self._sc_ser_layout.setSpacing(3); self._sc_ser_layout.setContentsMargins(0,0,0,0)
        self._sc_ser_layout.addStretch()
        sc_ser_scroll.setWidget(self._sc_ser_container)
        sc_ser_vl.addWidget(sc_ser_scroll, 1)
        sc_body.addWidget(sc_ser_frame)
        sl.addLayout(sc_body, 1)
        self.tabs.addTab(self.tab_summary, self.T["tab_summary"])
        self.tabs.addTab(self.tab_scatter, "Nokta Grafik")
        # Parametre/ölçek değişince yeniden çiz
        for btn in self._sc_param_btns.values():
            btn.toggled.connect(lambda c: self._on_scatter_param() if c else None)
        for btn in self._sc_scale_btns.values():
            btn.toggled.connect(lambda c: self._plot_scatter() if c else None)

        # Log-Probit
        self.tab_lp=QWidget()
        lpl=QVBoxLayout(self.tab_lp); lpl.setContentsMargins(4,4,4,4); lpl.setSpacing(4)
        lp_ctrl=QHBoxLayout()
        self.chk_lp_avg=QCheckBox(self.T["lp_avg_only"])
        self.chk_lp_avg.stateChanged.connect(
            lambda: self._plot_lp() if self.all_series else None)
        lp_ctrl.addWidget(self.chk_lp_avg); lp_ctrl.addStretch()
        lpl.addLayout(lp_ctrl)
        self.lp_canvas=FigureCanvas(Figure(figsize=(9,5), facecolor=BG))
        lpl.addWidget(self.lp_canvas,1)
        # Bilgi etiketi - grafiğin altında
        self.lp_info_lbl=QLabel("")
        self.lp_info_lbl.setVisible(False)
        self.lp_info_lbl.setTextFormat(Qt.TextFormat.RichText)
        self.lp_info_lbl.setWordWrap(True)
        self.lp_info_lbl.setStyleSheet(
            f"background:#111827;color:#90c0e0;font-size:11px;"
            f"padding:6px 10px;border-top:1px solid #1a2a40;border-radius:0 0 4px 4px;")
        lpl.addWidget(self.lp_info_lbl)
        self.tabs.addTab(self.tab_lp, self.T["tab_plot"])

        # Dağılım — tüm sekme scroll içinde
        self.tab_dist=QWidget()
        dist_outer=QVBoxLayout(self.tab_dist)
        dist_outer.setContentsMargins(0,0,0,0); dist_outer.setSpacing(0)
        dist_scroll=QScrollArea(); dist_scroll.setWidgetResizable(True)
        dist_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        dist_scroll.setStyleSheet("border:none;")
        dist_inner=QWidget()
        dl=QVBoxLayout(dist_inner); dl.setContentsMargins(4,4,4,4); dl.setSpacing(4)
        lim_panel=QFrame(); lim_panel.setFixedHeight(36)
        lim_panel.setStyleSheet(f"background:{BG3};border-radius:4px;border:1px solid #2a4060;")
        limhl=QHBoxLayout(lim_panel); limhl.setSpacing(8); limhl.setContentsMargins(8,2,8,2)
        limhl.addWidget(QLabel(self.T["limit_label"]))
        self.lim_grp=QButtonGroup(self)
        for val,txt in [("ema",self.T["lim_ema"]),("fda",self.T["lim_fda"]),
                         ("usp",self.T["lim_usp"]),("custom",self.T["lim_custom"])]:
            rb=QRadioButton(txt); rb.setProperty("lim_val",val)
            self.lim_grp.addButton(rb)
            if val=="ema": rb.setChecked(True)
            rb.toggled.connect(lambda c,v=val:
                (setattr(self,'limit_type',v) or self._plot_dist())
                if c and self.all_series else None)
            limhl.addWidget(rb)
        limhl.addWidget(QLabel(self.T["lim_pct"]))
        self.e_lim_pct=QLineEdit("20"); self.e_lim_pct.setFixedWidth(44)
        self.e_lim_pct.editingFinished.connect(
            lambda: self._plot_dist() if self.all_series else None)
        limhl.addWidget(self.e_lim_pct); limhl.addStretch()
        dl.addWidget(lim_panel)
        self.dist_canvas=FigureCanvas(Figure(figsize=(9,5), facecolor=BG))
        self.dist_canvas.setMinimumHeight(420)
        dl.addWidget(self.dist_canvas)
        # Uyarı scroll
        self.warn_scroll=QScrollArea(); self.warn_scroll.setWidgetResizable(True)
        self.warn_scroll.setMaximumHeight(160); self.warn_scroll.setVisible(False)
        self.warn_scroll.setStyleSheet("background:#1a0a0a;border:none;")
        self.warn_container=QWidget(); self.warn_container.setStyleSheet("background:#1a0a0a;")
        self.warn_layout=QVBoxLayout(self.warn_container)
        self.warn_layout.setContentsMargins(6,4,6,4); self.warn_layout.setSpacing(3)
        self.warn_scroll.setWidget(self.warn_container)
        dl.addWidget(self.warn_scroll)
        dist_inner.setLayout(dl)
        dist_scroll.setWidget(dist_inner)
        dist_outer.addWidget(dist_scroll)
        self.tabs.addTab(self.tab_dist, self.T["tab_dist"])

        # Özet
        self.tab_summary=QScrollArea(); self.tab_summary.setWidgetResizable(True)
        self.summary_widget=QWidget()
        self.summary_layout=QVBoxLayout(self.summary_widget)
        self.summary_layout.setSpacing(6); self.summary_layout.setContentsMargins(8,8,8,8)
        self.summary_layout.addStretch()
        self.tab_summary.setWidget(self.summary_widget)

        # Karşılaştırma
        self.tab_compare=QScrollArea(); self.tab_compare.setWidgetResizable(True)
        self.compare_widget=QWidget()
        self.compare_layout=QVBoxLayout(self.compare_widget)
        self.compare_layout.setSpacing(6); self.compare_layout.setContentsMargins(8,8,8,8)
        self.compare_layout.addStretch()
        self.tab_compare.setWidget(self.compare_widget)
        self.tabs.addTab(self.tab_compare, self.T["tab_compare"])

    # ── Yardımcı ──────────────────────────────────────────────────────────────
    def _refresh_cutoffs(self):
        for i in reversed(range(self.cbox_layout.count())):
            w=self.cbox_layout.itemAt(i).widget()
            if w: w.setParent(None)
        flow=int(self.flow_combo.currentText()); co=NGI_CUTOFFS[flow]
        lbl=QLabel(f"{flow}L:")
        lbl.setStyleSheet(f"color:{GOLD};font-size:9px;font-weight:bold;"
            "background:transparent;border:none;")
        self.cbox_layout.addWidget(lbl)
        for s in ["S1","S2","S3","S4","S5","S6","S7","MOC"]:
            if co.get(s,999)<900:
                b=QFrame()
                b.setStyleSheet("background:#1a2540;border-radius:2px;"
                    "border:0.5px solid #2a4060;")
                bl=QVBoxLayout(b); bl.setContentsMargins(3,1,3,1); bl.setSpacing(0)
                sl2=QLabel(s); sl2.setStyleSheet(
                    "color:#7ab0d0;font-size:8px;font-weight:bold;"
                    "background:transparent;border:none;")
                sl2.setAlignment(Qt.AlignmentFlag.AlignCenter)
                vl2=QLabel(f"{co[s]:.2f}"); vl2.setStyleSheet(
                    "color:#e0f0ff;font-size:8px;background:transparent;border:none;")
                vl2.setAlignment(Qt.AlignmentFlag.AlignCenter)
                bl.addWidget(sl2); bl.addWidget(vl2)
                self.cbox_layout.addWidget(b)
        self.cbox_layout.addStretch()

    def _on_flow(self, v):
        self.flow=int(v); self._refresh_cutoffs()

    def _add_series(self):
        idx=len(self.series_panels)+1
        color=CP[(idx-1)%len(CP)]
        p=SeriesPanel(idx, color, self.T, self.series_container)
        p.del_btn.clicked.connect(lambda: self._del_specific(p))
        p.vis_check.stateChanged.connect(lambda: self._replot_if_done())
        cnt=self.series_layout.count()
        self.series_layout.insertWidget(cnt-1, p)
        self.series_panels.append(p)

    def _del_series(self):
        if len(self.series_panels)<=1: return
        p=self.series_panels.pop()
        p.setParent(None)

    def _del_specific(self, p):
        if len(self.series_panels)<=1: return
        self.series_panels.remove(p)
        p.setParent(None)

    def _replot_if_done(self):
        if not self.all_series: return
        cur=self.tabs.currentIndex()
        if cur==1: self._plot_lp()
        elif cur==2: self._plot_dist()
        elif cur==4: self._show_compare()

    def _select_all_series(self):
        for p in self.series_panels: p.vis_check.setChecked(True)

    def _deselect_all_series(self):
        for p in self.series_panels: p.vis_check.setChecked(False)

    def _visible_series(self):
        return [sd for sd,p in zip(self.all_series, self.series_panels)
                if p.vis_check.isChecked()]

    def _toggle_lang(self):
        self.lang="EN" if self.lang=="TR" else "TR"
        self.T=L[self.lang]; T=self.T
        self.setWindowTitle(T["title"])
        self.lbl_title.setText(T["title"])
        self.lbl_sub.setText(T["subtitle"])
        self.btn_lang.setText(T["lang_btn"])
        self.btn_calc.setText(T["calculate"])
        self.btn_clr.setText(T["clear"])
        self.btn_pdf.setText(T["export_pdf"])
        self.btn_csv.setText(T["load_csv"])
        self.btn_add.setText("+ "+T["add_series"])
        self.chk_lp_avg.setText(T["lp_avg_only"])
        self.chk_delivered_tp.setText(T["delivered_tp"])
        try:
            self.btn_sel_all.setText("✓ "+("Tümü" if self.lang=="TR" else "All"))
            self.btn_sel_none.setText("✕ "+("Gizle" if self.lang=="TR" else "Hide"))
        except: pass
        self._refresh_cutoffs()
        for p in self.series_panels: p.update_lang(T)
        for i,k in enumerate(["tab_summary","tab_scatter","tab_plot","tab_dist","tab_compare"]):
            self.tabs.setTabText(i, T[k])
        self._update_series_count()
        if self.all_series:
            cur=self.tabs.currentIndex()
            if cur==0: self._show_summary()
            elif cur==1: self._plot_scatter()
            elif cur==2: self._plot_lp()
            elif cur==3: self._plot_dist()
            elif cur==4: self._show_compare()

    def _clear(self):
        for p in self.series_panels:
            for ri in range(RUNS_PER_SER):
                for s in DISP_STAGES:
                    p.entries[ri][s].setText("0.000")
            if p.ref_check: p.ref_check.setChecked(False)
        self.all_series=[]
        self._clear_results()

    def _clear_results(self):
        for lay in [self.summary_layout, self.compare_layout]:
            for i in reversed(range(lay.count()-1)):
                w=lay.itemAt(i).widget()
                if w: w.setParent(None)
        for canvas in [self.lp_canvas, self.dist_canvas, self.scatter_canvas]:
            canvas.figure.clear(); canvas.draw()
        for i in reversed(range(self.warn_layout.count())):
            w=self.warn_layout.itemAt(i).widget()
            if w: w.setParent(None)
        self.warn_scroll.setVisible(False)

    def _get_flow_lo_hi(self):
        flow=int(self.flow_combo.currentText())
        try: lo=float(self.e_lo.text())
        except: lo=15
        try: hi=float(self.e_hi.text())
        except: hi=85
        return flow,lo,hi

    # ── Hesaplama ─────────────────────────────────────────────────────────────
    def _calculate(self):
        flow,lo,hi=self._get_flow_lo_hi()
        si=[{"name":p.name_edit.text(),"color":p.color,
             "masses_list":[p.get_masses(ri) for ri in range(RUNS_PER_SER)],
             "is_ref":(p.ref_check.isChecked() if p.ref_check else False)
                       and p==self.series_panels[0]}
            for p in self.series_panels]
        self.status_lbl.setText(self.T["status_calc"])
        self.btn_calc.setEnabled(False)
        self.calc_thread=CalcThread(si,flow,lo,hi,self.chk_delivered_tp.isChecked())
        self.calc_thread.done.connect(self._on_calc_done)
        self.calc_thread.error.connect(self._on_calc_error)
        self.calc_thread.finished.connect(lambda: self._thread_cleanup())
        self.calc_thread.start()
        # Thread referansını listede tut (GC'den koru)
        if not hasattr(self,'_threads'): self._threads=[]
        self._threads.append(self.calc_thread)

    def _on_calc_done(self, results):
        self.all_series=results; self.btn_calc.setEnabled(True)
        self._clear_results()
        self._rebuild_scatter_series_list()
        cur=self.tabs.currentIndex()
        if cur==0: self._show_summary()
        elif cur==1: self._plot_scatter()
        elif cur==2: self._plot_lp()
        elif cur==3: self._plot_dist()
        elif cur==4: self._show_compare()
        self.status_lbl.setText(self.T["status_done"])

    def _on_calc_error(self, msg):
        self.btn_calc.setEnabled(True)
        QMessageBox.critical(self,"Hesaplama Hatasi",msg[:500])

    def _thread_cleanup(self):
        """Tamamlanan thread'leri listeden temizle"""
        if hasattr(self,'_threads'):
            self._threads=[t for t in self._threads if t.isRunning()]

    def _on_tab_change(self, idx):
        if not self.all_series: return
        if idx==0: self._show_summary()
        elif idx==1: self._plot_scatter()
        elif idx==2: self._plot_lp()
        elif idx==3: self._plot_dist()
        elif idx==4: self._show_compare()

    # ── Sonuçlar ──────────────────────────────────────────────────────────────
    # ── Nokta Grafik metodları ────────────────────────────────────────────────
    def _rebuild_scatter_series_list(self):
        """Seri listesini yeniden oluştur"""
        for i in reversed(range(self._sc_ser_layout.count()-1)):
            w = self._sc_ser_layout.itemAt(i).widget()
            if w: w.setParent(None)
        self._sc_ser_checks = {}
        for sd in self.all_series:
            row = QWidget()
            row.setStyleSheet("background:transparent;")
            rl = QHBoxLayout(row); rl.setContentsMargins(0,1,0,1); rl.setSpacing(4)
            bar = QFrame(); bar.setFixedSize(3,16)
            bar.setStyleSheet(f"background:{sd['color']};border-radius:1px;")
            rl.addWidget(bar)
            chk = QCheckBox(sd["name"])
            chk.setChecked(True)
            chk.setStyleSheet(f"""
                QCheckBox{{color:#c0d8f0;font-size:11px;spacing:4px;}}
                QCheckBox::indicator{{width:13px;height:13px;
                    border:1.5px solid {sd['color']};border-radius:2px;background:transparent;}}
                QCheckBox::indicator:checked{{background:{sd['color']};}}
            """)
            chk.stateChanged.connect(lambda _: self._plot_scatter())
            rl.addWidget(chk, 1)
            self._sc_ser_layout.insertWidget(
                self._sc_ser_layout.count()-1, row)
            self._sc_ser_checks[sd["name"]] = chk

    def _scatter_sel_all(self):
        for chk in self._sc_ser_checks.values():
            chk.setChecked(True)

    def _scatter_sel_none(self):
        for chk in self._sc_ser_checks.values():
            chk.setChecked(False)

    def _on_scatter_param(self):
        """Parametre butonu değişince alt ölçek görünürlüğünü ayarla"""
        param = self._get_scatter_param()
        show_scale = param in ("fpd","fpf")
        self.sc_row2_widget.setVisible(show_scale)
        self._plot_scatter()

    

    def _get_scatter_param(self):
        for key, btn in self._sc_param_btns.items():
            if btn.isChecked(): return key
        return "fpd"

    def _get_scatter_scale(self):
        for key, btn in self._sc_scale_btns.items():
            if btn.isChecked(): return key
        return "5"

    def _get_scatter_key(self):
        """Gerçek parametre anahtarını döndür"""
        param = self._get_scatter_param()
        scale = self._get_scatter_scale()
        if param == "fpd":
            return {"5":"fpd","3":"fpd3","15":"fpd15"}[scale]
        elif param == "fpf":
            return {"5":"fpf","3":"fpf3","15":"fpf15"}[scale]
        else:
            return param  # mmad, gsd, d10, d50, d90

    def _get_scatter_label(self):
        param = self._get_scatter_param()
        scale = self._get_scatter_scale()
        if param == "fpd":
            d = {"5":"5","3":"3","15":"1,5"}[scale]
            return f"FPD ≤{d} µm [mg]"
        elif param == "fpf":
            d = {"5":"5","3":"3","15":"1,5"}[scale]
            return f"FPF ≤{d} µm [%]"
        elif param == "mmad": return "MMAD [µm]"
        elif param == "gsd":  return "GSD"
        elif param == "d10":  return "D10 [µm]"
        elif param == "d50":  return "D50 [µm]"
        elif param == "d90":  return "D90 [µm]"
        return param

    def _plot_scatter(self):
        if not self.all_series: return
        plt.close("all")
        fig=self.scatter_canvas.figure; fig.clear()
        fig.patch.set_facecolor(BG)
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")

        key=self._get_scatter_key()
        lbl=self._get_scatter_label()
        ds=self.T["dec_sep"]

        lm={"ema":20,"fda":15,"usp":25}
        try: pct=lm.get(self.limit_type) or float(self.e_lim_pct.text())
        except: pct=20

        ref_val=None
        ref_sd=next((sd for sd in self.all_series if sd["is_ref"]),None)
        if ref_sd and ref_sd["avg"] and key in ref_sd["avg"]["params"]:
            ref_val=ref_sd["avg"]["params"][key][0]

        visible=[sd for sd in self.all_series
                 if sd["name"] in getattr(self,"_sc_ser_checks",{})
                 and self._sc_ser_checks[sd["name"]].isChecked()]

        if not visible:
            ax.text(0.5,0.5,"Seri seçilmedi",transform=ax.transAxes,
                ha="center",va="center",color=TXT2,fontsize=13)
            fig.tight_layout(); self.scatter_canvas.draw(); return

        n=len(visible)
        # Limit çizgileri
        if ref_val is not None and ref_val>0:
            ax.axhline(ref_val,color="#4a9ae8",lw=2.0,zorder=2,
                label=f"Ref = {fmt_num(ref_val,3,ds)}")
            up=ref_val*(1+pct/100); lo=ref_val*(1-pct/100)
            ax.axhline(up,color="#FF4444",lw=1.5,ls="--",zorder=2,
                label=f"+{pct:.0f}%")
            ax.axhline(lo,color="#FF4444",lw=1.5,ls="--",zorder=2,
                label=f"–{pct:.0f}%")
            ax.fill_between([-0.5,n-0.5],lo,up,color="#FF4444",alpha=0.05,zorder=1)

        # Değerleri topla (Y limiti hesabı için)
        all_vals=[sd["avg"]["params"][key][0]
                  for sd in visible
                  if sd["avg"] and key in sd["avg"]["params"]]
        if not all_vals:
            fig.tight_layout(); self.scatter_canvas.draw(); return

        y_min=min(all_vals); y_max=max(all_vals)
        if ref_val:
            y_min=min(y_min, ref_val*(1-pct/100)*0.9)
            y_max=max(y_max, ref_val*(1+pct/100)*1.1)
        y_range=max(y_max-y_min, y_max*0.1, 1e-6)
        y_pad=y_range*0.25

        # Noktaları çiz
        for xi,sd in enumerate(visible):
            if not sd["avg"] or key not in sd["avg"]["params"]: continue
            val,val_sd=sd["avg"]["params"][key][0],sd["avg"]["params"][key][1]
            is_ref=sd["is_ref"]
            marker="*" if is_ref else "o"
            ms=18 if is_ref else 10
            ax.plot(xi,val,marker=marker,color=sd["color"],
                markersize=ms,markeredgecolor="white",markeredgewidth=0.8,
                zorder=5 if is_ref else 4)
            # SD çubuğu
            if val_sd>0:
                ax.errorbar(xi,val,yerr=val_sd,fmt="none",
                    color=sd["color"],capsize=5,lw=1.5,alpha=0.7,zorder=3)
            # Değer etiketi (noktanın üstünde)
            ax.annotate(fmt_num(val,3,ds),(xi,val),
                xytext=(0,12),textcoords="offset points",
                ha="center",va="bottom",fontsize=10,
                color=sd["color"],fontweight="bold" if is_ref else "normal")
            # Seri adı etiketi (noktanın altında, kısa)
            short_name=sd["name"][:12]+"…" if len(sd["name"])>12 else sd["name"]
            ax.annotate(short_name,(xi,val),
                xytext=(0,-16),textcoords="offset points",
                ha="center",va="top",fontsize=9,
                color="#8090b0",rotation=0)
            # Limit dışı kırmızı çerçeve
            if ref_val and ref_val>0:
                if val>ref_val*(1+pct/100) or val<ref_val*(1-pct/100):
                    ax.plot(xi,val,"o",color="none",
                        markeredgecolor="#FF4444",markeredgewidth=2.5,
                        markersize=20,zorder=6)

        # Eksen düzenlemeleri
        ax.set_xticks(list(range(n)))
        ax.set_xticklabels([""]*n)  # Seri adları annotate ile yazıldı
        ax.set_xlim(-0.6,n-0.4)
        ax.set_ylim(y_min-y_pad, y_max+y_pad)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6,prune="both"))
        ax.tick_params(colors=TXT2,labelsize=11)
        for sp in ax.spines.values(): sp.set_color("#2a4060")
        ax.set_ylabel(lbl,color=TXT2,fontsize=12)
        title=f"Nokta Grafik — {lbl}"
        if ref_val: title+=f"  |  Limit ±{pct:.0f}%"
        ax.set_title(title,color=GOLD,fontsize=13,fontweight="bold")
        ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
        hdls,_=ax.get_legend_handles_labels()
        if hdls:
            ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0",
                framealpha=0.85,loc="upper right")
        fig.tight_layout(); self.scatter_canvas.draw()
    # ── Log-Probit ────────────────────────────────────────────────────────────
    def _plot_lp(self):
        plt.close("all")
        fig=self.lp_canvas.figure; fig.clear(); fig.patch.set_facecolor(BG)
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        flow=int(self.flow_combo.currentText())
        avg_only=self.chk_lp_avg.isChecked() or len(self.all_series)>=4
        if avg_only and not self.chk_lp_avg.isChecked():
            self.chk_lp_avg.blockSignals(True); self.chk_lp_avg.setChecked(True)
            self.chk_lp_avg.blockSignals(False)
        for sd in self._visible_series():
            col=sd["color"]
            vr=[r for r in sd["runs"] if "error" not in r]
            if not vr: continue
            lw=2.5 if sd["is_ref"] else 1.5
            if avg_only:
                ml=min(len(r["x_reg"]) for r in vr)
                ax_=sum(r["x_reg"][:ml] for r in vr)/len(vr)
                ay_=sum(r["y_reg"][:ml] for r in vr)/len(vr)
                ba=sum(r["b"] for r in vr)/len(vr)
                aa=sum(r["a"] for r in vr)/len(vr)
                ax.plot(ax_,ay_,"o",color=col,ms=6,zorder=4)
                xr=np.linspace(min(ax_)-0.15,max(ax_)+0.15,60)
                ax.plot(xr,aa+ba*xr,"-",color=col,lw=lw,label=sd["name"])
                if sd["avg"] and "mmad" in sd["avg"]["params"]:
                    mv=sd["avg"]["params"]["mmad"][0]
                    if mv>0: ax.axvline(math.log10(mv),color=col,lw=0.8,ls=":",alpha=0.7)
            else:
                for run in vr:
                    ax.plot(run["x_reg"],run["y_reg"],"o",color=col,ms=4,alpha=0.7)
                    xr=np.linspace(min(run["x_reg"])-0.1,max(run["x_reg"])+0.1,50)
                    ax.plot(xr,run["a"]+run["b"]*xr,"-",color=col,lw=lw,alpha=0.7,
                        label=f"{sd['name']} R{run['run_no']}")
        ax.set_xlabel("log₁₀(D50, µm)",color=TXT2,fontsize=13)
        ax.set_ylabel("Probit z",color=TXT2,fontsize=13)
        ax.set_title(f"Log-Probit  [{flow} L/min]",color=GOLD,fontsize=14,fontweight="bold")
        ax.tick_params(colors=TXT2,labelsize=11)
        for sp in ax.spines.values(): sp.set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        hdls,_=ax.get_legend_handles_labels()
        if hdls:
            n=len(hdls)
            ax.legend(fontsize=max(9,12-max(0,n-6)),facecolor="#0e1525",
                labelcolor="#d0e0f0",ncol=2 if n>6 else 1,framealpha=0.85)
        fig.tight_layout(); self.lp_canvas.draw()
        # Bilgi tablosu altta
        vis=[sd for sd in self._visible_series()
             if sd.get("avg") and "mmad" in sd["avg"]["params"]]
        if vis:
            ds=self.T["dec_sep"]; info_lines=[]
            for sd in vis:
                pr=sd["avg"]["params"]
                info_lines.append(
                    f"<span style='color:{sd['color']};font-weight:bold;'>{sd['name']}</span>"
                    f"&nbsp;&nbsp;MMAD={fmt_num(pr['mmad'][0],3,ds)}"
                    f"&nbsp;GSD={fmt_num(pr['gsd'][0],3,ds)}"
                    f"&nbsp;D10={fmt_num(pr.get('d10',(0,))[0],3,ds)}"
                    f"&nbsp;D50={fmt_num(pr.get('d50',(0,))[0],3,ds)}"
                    f"&nbsp;D90={fmt_num(pr.get('d90',(0,))[0],3,ds)}"
                    f"&nbsp;Slope={fmt_num(pr['slope'][0],3,ds)}"
                    f"&nbsp;R²={fmt_num(pr['r2'][0],4,ds)}")
            try:
                self.lp_info_lbl.setText("<br>".join(info_lines))
                self.lp_info_lbl.setVisible(True)
            except: pass
        else:
            try: self.lp_info_lbl.setVisible(False)
            except: pass

    # ── APSD Dağılım ──────────────────────────────────────────────────────────
    def _plot_dist(self):
        plt.close("all")
        flow=int(self.flow_combo.currentText())
        co=NGI_CUTOFFS[flow]
        vis_all=["Throat"]+[s for s in GRAPH_STAGES if co.get(s,999)<900]
        x_all=np.arange(len(vis_all))
        lm={"ema":20,"fda":15,"usp":25}
        try: pct=lm.get(self.limit_type) or float(self.e_lim_pct.text())
        except: pct=20
        fig=self.dist_canvas.figure; fig.clear(); fig.patch.set_facecolor(BG)
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        ref_masses=None; warnings=[]
        for sd in self._visible_series():
            if not sd["avg"]: continue
            ms=[sd["avg"]["avg_masses"].get(s,0) for s in vis_all]
            vr=[r for r in sd["runs"] if "error" not in r]
            sds=[]
            for s in vis_all:
                vals=[r["masses"].get(s,0) for r in vr]
                sds.append(float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0)
            lw=3 if sd["is_ref"] else 1.8
            ax.plot(x_all,ms,color=sd["color"],lw=lw,marker="o",
                markersize=10 if sd["is_ref"] else 5,label=sd["name"],zorder=4)
            ax.fill_between(x_all,[m-s for m,s in zip(ms,sds)],
                [m+s for m,s in zip(ms,sds)],color=sd["color"],alpha=0.12,zorder=2)
            ax.errorbar(x_all,ms,yerr=sds,fmt="none",
                color=sd["color"],capsize=4,lw=1.5,alpha=0.5,zorder=3)
            if sd["is_ref"]: ref_masses=sd["avg"]["avg_masses"]
        if ref_masses:
            rv=[ref_masses.get(s,0) for s in vis_all]
            upper=[v*(1+pct/100) for v in rv]; lower=[v*(1-pct/100) for v in rv]
            ax.plot(x_all,upper,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"+{pct:.0f}%")
            ax.plot(x_all,lower,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"–{pct:.0f}%")
            ax.fill_between(x_all,lower,upper,color="#FF6060",alpha=0.05,zorder=1)
            for sd in self._visible_series():
                if sd["is_ref"] or not sd["avg"]: continue
                mt=[sd["avg"]["avg_masses"].get(s,0) for s in vis_all]
                for s,tv,lo2,hi2 in zip(vis_all,mt,lower,upper):
                    if tv<lo2 or tv>hi2:
                        yon=("yuksek" if tv>hi2 else "dusuk") if self.lang=="TR" else ("high" if tv>hi2 else "low")
                        warnings.append((f"{sd['name']} - {s}: {fmt_num(tv,4,self.T['dec_sep'])} ({yon})",False,False))
                f2=calc_f2(ref_masses,sd["avg"]["avg_masses"],co)
                if f2 is not None:
                    pf2=self.T["f2_pass"] if f2>=50 else self.T["f2_fail"]
                    warnings.insert(0,(f"{self.T['f2_label']} {sd['name']}: f2={fmt_num(f2,1,self.T['dec_sep'])} ({pf2})",f2>=50,True))
        ax.set_xticks(list(x_all))
        ax.set_xticklabels(vis_all,rotation=20,ha="right",fontsize=12)
        ax.tick_params(colors=TXT2,labelsize=11)
        for sp in ax.spines.values(): sp.set_color("#2a4060")
        ylbl="Ort. Kutle (mg/atis)" if self.lang=="TR" else "Mean Mass (mg/actuation)"
        ax.set_xlabel(self.T["stage"],color=TXT2,fontsize=13)
        ax.set_ylabel(ylbl,color=TXT2,fontsize=13)
        ttl=f"APSD [{flow} L/min] Ort±SD"
        if ref_masses: ttl+=f"  |  Limit ±{pct:.0f}%"
        ax.set_title(ttl,color=GOLD,fontsize=13,fontweight="bold")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        hdls,_=ax.get_legend_handles_labels()
        if hdls:
            n=len(hdls)
            ax.legend(fontsize=max(9,12-max(0,n-6)),facecolor="#0e1525",
                labelcolor="#d0e0f0",ncol=2 if n>6 else 1,framealpha=0.85,loc="upper right")
        fig.tight_layout(); self.dist_canvas.draw()
        for i in reversed(range(self.warn_layout.count())):
            w=self.warn_layout.itemAt(i).widget()
            if w: w.setParent(None)
        if warnings:
            self.warn_scroll.setVisible(True)
            for item in warnings:
                wt,is_pass,is_f2=item
                bg=("#0a2a0a" if is_pass else "#2a0a0a") if is_f2 else "#2a0a0a"
                tc=("#90ee90" if is_pass else "#FF6060") if is_f2 else "#FFB0B0"
                wl=QLabel(f"  {wt}")
                wl.setStyleSheet(f"color:{tc};font-weight:bold;font-size:12px;"
                    f"background:{bg};border-radius:3px;padding:3px 6px;")
                self.warn_layout.addWidget(wl)
        else: self.warn_scroll.setVisible(False)

    # ── Özet / Sonuçlar ───────────────────────────────────────────────────────
    def _show_summary(self):
        for i in reversed(range(self.summary_layout.count()-1)):
            w=self.summary_layout.itemAt(i).widget()
            if w: w.setParent(None)
        T=self.T; ds=T["dec_sep"]
        try: rsd_lim=float(self.e_rsd.text())
        except: rsd_lim=5.0
        params_list=[
            ("metered",T["metered"]),("delivered",T["delivered"]),
            ("fpd",T["fp_dose"]+" (5µm)"),("fpf",T["fp_frac"]+" (5µm)"),
            ("fpd3","FPD (3µm) [mg]"),("fpf3","FPF (3µm) [%]"),
            ("fpd15","FPD (1.5µm) [mg]"),("fpf15","FPF (1.5µm) [%]"),
            ("mmad","MMAD (µm)"),("gsd","GSD"),
            ("d10","D10 (µm)"),("d50","D50 (µm)"),("d90","D90 (µm)"),
            ("slope",T["slope_lbl"]),("intercept",T["int_lbl"]),("r2",T["r2_lbl"])]
        for sd in self.all_series:
            rt=f"  [{T['ref_label']}]" if sd["is_ref"] else ""
            sh=QLabel(f"▌ {sd['name']}{rt}")
            sh.setStyleSheet(f"color:{sd['color']};font-size:14px;font-weight:bold;"
                f"background:{BG2};border-left:3px solid {sd['color']};padding:4px 8px;")
            self.summary_layout.insertWidget(self.summary_layout.count()-1,sh)
            vr=[r for r in sd["runs"] if "error" not in r]
            if not vr: continue
            n=len(vr)
            tf=QFrame(); tf.setStyleSheet("background:#111827;border-radius:4px;border:1px solid #1a2a40;")
            tfl=QVBoxLayout(tf); tfl.setContentsMargins(2,2,2,2); tfl.setSpacing(1)
            ws=[150]+[90]*n+[90,90,76,62]
            hdrs=[T["param"]]+[f"Run {r['run_no']}" for r in vr]+[T["mean"],T["sd"],T["rsd"],T["accept"]]
            hf=QFrame(); hf.setStyleSheet("background:#1F4E79;border-radius:2px;")
            hfl=QHBoxLayout(hf); hfl.setContentsMargins(2,2,2,2); hfl.setSpacing(0)
            for h,w in zip(hdrs,ws):
                l=QLabel(h); l.setFixedWidth(w); l.setAlignment(Qt.AlignmentFlag.AlignCenter)
                l.setStyleSheet("color:white;font-weight:bold;font-size:11px;background:transparent;")
                hfl.addWidget(l)
            hfl.addStretch(); tfl.addWidget(hf)
            for i2,(key,lbl) in enumerate(params_list):
                vals=[r.get(key) for r in vr if key in r and r.get(key) is not None]
                if not vals: continue
                mv=float(np.mean(vals)); sv=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
                rv2=sv/mv*100 if mv else 0.0; pf=rv2<=rsd_lim
                ik=key in("fpd","fpf","fpd3","fpf3","fpd15","fpf15","mmad","gsd","d10","d50","d90")
                bg="#1a1a2a" if i2%2==0 else "transparent"
                dr=QFrame(); dr.setStyleSheet(f"background:{bg};")
                dfl=QHBoxLayout(dr); dfl.setContentsMargins(2,1,2,1); dfl.setSpacing(0)
                l0=QLabel(lbl); l0.setFixedWidth(ws[0])
                l0.setStyleSheet(f"color:{'#FFC600' if ik else '#c0d0e0'};"
                    f"font-weight:{'bold' if ik else 'normal'};font-size:11px;background:transparent;")
                dfl.addWidget(l0)
                for r2 in vr:
                    v=r2.get(key) or 0
                    l2=QLabel(fmt_num(v,4,ds)); l2.setFixedWidth(90)
                    l2.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    l2.setStyleSheet("color:#d0e8ff;font-size:11px;background:transparent;")
                    dfl.addWidget(l2)
                for val,w in [(fmt_num(mv,4,ds),90),(fmt_num(sv,4,ds),90),(fmt_num(rv2,2,ds),76)]:
                    l3=QLabel(val); l3.setFixedWidth(w); l3.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    l3.setStyleSheet("color:#d0e8ff;font-size:11px;background:transparent;")
                    dfl.addWidget(l3)
                pl=QLabel("OK" if pf else "FAIL"); pl.setFixedWidth(62)
                pl.setAlignment(Qt.AlignmentFlag.AlignCenter)
                pl.setStyleSheet(f"color:{'#90ee90' if pf else '#ff6060'};"
                    "font-weight:bold;font-size:11px;background:transparent;")
                dfl.addWidget(pl); dfl.addStretch(); tfl.addWidget(dr)
            self.summary_layout.insertWidget(self.summary_layout.count()-1,tf)
            dv=[r.get("delivered",0) for r in vr]
            if dv:
                dm=float(np.mean(dv)); ds2=float(np.std(dv,ddof=1)) if len(dv)>1 else 0.0
                dr2=ds2/dm*100 if dm else 0.0
                dl=QLabel(f"  {T['ddu_label']}: {T['mean']}={fmt_num(dm,4,ds)}mg  "
                    f"{T['sd']}={fmt_num(ds2,4,ds)}  {T['rsd']}={fmt_num(dr2,2,ds)}%")
                dl.setStyleSheet("color:#90ee90;font-size:12px;background:#1a2a1a;"
                    "border-radius:4px;padding:4px 8px;")
                self.summary_layout.insertWidget(self.summary_layout.count()-1,dl)

    # ── Karşılaştırma ─────────────────────────────────────────────────────────
    def _show_compare(self):
        for i in reversed(range(self.compare_layout.count()-1)):
            w=self.compare_layout.itemAt(i).widget()
            if w: w.setParent(None)
        T=self.T; ds=T["dec_sep"]
        flow=int(self.flow_combo.currentText()); co=NGI_CUTOFFS[flow]
        vis_ser=self._visible_series()
        if len(vis_ser)>=2:
            plt.close("all")
            fig2=Figure(figsize=(9,3.2),facecolor=BG)
            canvas2=FigureCanvas(fig2); canvas2.setMaximumHeight(260)
            ax1=fig2.add_subplot(121); ax2=fig2.add_subplot(122)
            ax1.set_facecolor("#0e1525"); ax2.set_facecolor("#0e1525")
            vsd=[sd for sd in vis_ser if sd["avg"]]
            names=[sd["name"] for sd in vsd]
            mmads=[sd["avg"]["params"].get("mmad",(0,))[0] for sd in vsd]
            gsds=[sd["avg"]["params"].get("gsd",(0,))[0] for sd in vsd]
            clrs=[sd["color"] for sd in vsd]
            xi=range(len(names))
            for i,(m,g,c) in enumerate(zip(mmads,gsds,clrs)):
                ax1.bar(i,m,color=c,alpha=0.85,width=0.6)
                ax2.bar(i,g,color=c,alpha=0.85,width=0.6)
            for ax,ttl in [(ax1,T["trend_mmad"]),(ax2,T["trend_gsd"])]:
                ax.set_xticks(list(xi))
                ax.set_xticklabels(names,rotation=20,ha="right",fontsize=9,color="#c0d8f0")
                ax.set_title(ttl,color=GOLD,fontsize=11,fontweight="bold")
                ax.tick_params(colors=TXT2)
                for sp in ax.spines.values(): sp.set_color("#2a4060")
                ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
            fig2.tight_layout()
            self.compare_layout.insertWidget(self.compare_layout.count()-1,canvas2)
        ref_sd=next((sd for sd in vis_ser if sd["is_ref"]),None)
        params_list=[
            ("mmad","MMAD (µm)"),("gsd","GSD"),
            ("d10","D10 (µm)"),("d50","D50 (µm)"),("d90","D90 (µm)"),
            ("fpd",T["fp_dose"]+" (5µm)"),("fpf",T["fp_frac"]+" (5µm)"),
            ("fpd3","FPD (3µm)"),("fpf3","FPF (3µm)"),
            ("fpd15","FPD (1.5µm)"),("fpf15","FPF (1.5µm)"),
            ("slope",T["slope_lbl"]),("intercept",T["int_lbl"]),("r2",T["r2_lbl"])]
        tf=QFrame(); tf.setStyleSheet("background:#111827;border-radius:4px;border:1px solid #1a2a40;")
        tfl=QVBoxLayout(tf); tfl.setContentsMargins(2,2,2,2); tfl.setSpacing(1)
        n_ser=len(vis_ser); ws=[140]+[104]*n_ser
        hdrs=[T["param"]]+[sd["name"]+(" ★" if sd["is_ref"] else "") for sd in vis_ser]
        hf=QFrame(); hf.setStyleSheet("background:#1F4E79;border-radius:2px;")
        hfl=QHBoxLayout(hf); hfl.setContentsMargins(2,2,2,2); hfl.setSpacing(0)
        for h,w in zip(hdrs,ws):
            l=QLabel(h); l.setFixedWidth(w); l.setAlignment(Qt.AlignmentFlag.AlignCenter)
            l.setStyleSheet("color:white;font-weight:bold;font-size:11px;background:transparent;")
            hfl.addWidget(l)
        hfl.addStretch(); tfl.addWidget(hf)
        for i2,(key,lbl) in enumerate(params_list):
            bg="#1a1a2a" if i2%2==0 else "transparent"
            dr=QFrame(); dr.setStyleSheet(f"background:{bg};")
            dfl=QHBoxLayout(dr); dfl.setContentsMargins(2,1,2,1); dfl.setSpacing(0)
            ik=key in("mmad","gsd","d10","d50","d90","fpd","fpf","fpd3","fpf3","fpd15","fpf15")
            l0=QLabel(lbl); l0.setFixedWidth(ws[0])
            l0.setStyleSheet(f"color:{'#FFC600' if ik else '#c0d0e0'};"
                f"font-weight:{'bold' if ik else 'normal'};font-size:11px;background:transparent;")
            dfl.addWidget(l0)
            rv=None
            if ref_sd and ref_sd["avg"] and key in ref_sd["avg"]["params"]:
                rv=ref_sd["avg"]["params"][key][0]
            for sd in vis_ser:
                if not sd["avg"] or key not in sd["avg"]["params"]:
                    l2=QLabel("-"); l2.setFixedWidth(104); l2.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    l2.setStyleSheet("color:#888;font-size:11px;background:transparent;")
                    dfl.addWidget(l2); continue
                val=sd["avg"]["params"][key][0]; txt=fmt_num(val,4,ds); clr="#d0e8ff"
                if rv and not sd["is_ref"] and rv>0:
                    diff=(val-rv)/rv*100
                    diff_str=fmt_num(abs(diff),1,ds); sign="+" if diff>=0 else "-"
                    txt+=f"\n({sign}{diff_str}%)"
                    clr="#90ee90" if abs(diff)<10 else "#FFB060" if abs(diff)<20 else "#FF6060"
                l2=QLabel(txt); l2.setFixedWidth(104); l2.setAlignment(Qt.AlignmentFlag.AlignCenter)
                l2.setStyleSheet(f"color:{clr};font-size:11px;background:transparent;")
                dfl.addWidget(l2)
            dfl.addStretch(); tfl.addWidget(dr)
        self.compare_layout.insertWidget(self.compare_layout.count()-1,tf)
        if ref_sd and ref_sd["avg"]:
            f2h=QLabel(f"  {T['f2_label']}")
            f2h.setStyleSheet(f"color:#FFD700;font-size:14px;font-weight:bold;"
                f"background:{BG2};border-left:3px solid #FFD700;padding:4px 8px;")
            self.compare_layout.insertWidget(self.compare_layout.count()-1,f2h)
            lm={"ema":20,"fda":15,"usp":25}
            try: pct=lm.get(self.limit_type) or float(self.e_lim_pct.text())
            except: pct=20
            for sd in vis_ser:
                if sd["is_ref"] or not sd["avg"]: continue
                f2=calc_f2(ref_sd["avg"]["avg_masses"],sd["avg"]["avg_masses"],co)
                if f2 is None: continue
                pf2=T["f2_pass"] if f2>=50 else T["f2_fail"]
                clr="#90ee90" if f2>=50 else "#FF6060"
                bg_f2="#0a2a0a" if f2>=50 else "#2a0a0a"
                ff=QLabel(f"  {sd['name']}  f2 = {fmt_num(f2,1,ds)}   {pf2}")
                ff.setStyleSheet(f"color:{clr};font-size:13px;font-weight:bold;"
                    f"background:{bg_f2};border-radius:4px;padding:6px 8px;")
                self.compare_layout.insertWidget(self.compare_layout.count()-1,ff)

    # ── CSV Yükle ─────────────────────────────────────────────────────────────
    def _load_csv(self):
        path,_=QFileDialog.getOpenFileName(self,"CSV Sec","",
            "CSV (*.csv);;All Files (*.*)")
        if not path: return
        try: series_dict,flow,warnings=parse_csv(path)
        except Exception as ex:
            QMessageBox.critical(self,"CSV Hatasi",str(ex)); return
        if series_dict is None:
            QMessageBox.critical(self,"CSV Hatasi",
                warnings[0] if warnings else "Okunamadi"); return
        for w in warnings:
            if w=="csv_err_flow": QMessageBox.warning(self,"",self.T["csv_err_flow"])
            elif w.startswith("csv_4runs__"):
                QMessageBox.warning(self,"",self.T["csv_4runs"].format(
                    s=w.replace("csv_4runs__","")))
        has_ref=any(v["ref"] for v in series_dict.values())
        if not has_ref:
            dlg=QDialog(self); dlg.setWindowTitle(self.T["csv_ref_ask"])
            dlg.setMinimumWidth(340); dlg.setStyleSheet(STYLE)
            vl2=QVBoxLayout(dlg); vl2.addWidget(QLabel(self.T["csv_ref_ask"]))
            lst=QListWidget()
            lst.setStyleSheet(f"background:{BG3};border:1px solid #2a4060;")
            for name in list(series_dict.keys())+[self.T["csv_ref_none"]]:
                lst.addItem(QListWidgetItem(name))
            lst.setCurrentRow(0); vl2.addWidget(lst)
            btns=QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel)
            btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
            vl2.addWidget(btns)
            if dlg.exec()==QDialog.DialogCode.Accepted:
                sel=lst.currentItem().text() if lst.currentItem() else ""
                if sel and sel!=self.T["csv_ref_none"]:
                    for k in series_dict: series_dict[k]["ref"]=(k==sel)
        if flow in NGI_CUTOFFS:
            self.flow_combo.setCurrentText(str(flow)); self._on_flow(str(flow))
        for p in self.series_panels: p.setParent(None)
        self.series_panels.clear()
        self.all_series=[]; self._clear_results()
        for si,(name,data) in enumerate(series_dict.items()):
            self._add_series()
            p=self.series_panels[-1]; p.name_edit.setText(name)
            if si==0 and p.ref_check: p.ref_check.setChecked(data["ref"])
            for ri,run_data in enumerate(data["runs"][:RUNS_PER_SER]):
                p.set_masses(ri,run_data["masses"])
        n_s=len(series_dict)
        n_r=sum(len(v["runs"]) for v in series_dict.values())
        self._update_series_count()
        self.status_lbl.setText(self.T["csv_loaded"].format(n=n_s,r=n_r))

    # ── PDF Rapor ─────────────────────────────────────────────────────────────
    def _export_pdf(self):
        if not self.all_series:
            QMessageBox.warning(self,"","Oncelikle hesaplama yapiniz."); return
        path,_=QFileDialog.getSaveFileName(self,"PDF Kaydet",
            f"NGI_{datetime.datetime.now().strftime('%Y%m%d_%H%M')}.pdf","PDF (*.pdf)")
        if not path: return
        meta={"product":self.e_product.text(),"batch":self.e_batch.text(),
              "operator":self.e_operator.text(),"date":self.e_date.text()}
        lm={"ema":20,"fda":15,"usp":25}
        try: pct=lm.get(self.limit_type) or float(self.e_lim_pct.text())
        except: pct=20
        try: rsd_lim=float(self.e_rsd.text())
        except: rsd_lim=5.0
        try:
            make_pdf_multi(path,self.all_series,meta,
                int(self.flow_combo.currentText()),
                self.T,pct,rsd_lim,lang=self.lang)
            import subprocess,platform
            try:
                if platform.system()=="Windows": os.startfile(path)
                elif platform.system()=="Darwin": subprocess.Popen(["open",path])
                else: subprocess.Popen(["xdg-open",path])
            except: pass
            QMessageBox.information(self,"",f"PDF kaydedildi:\n{path}")
        except Exception as ex:
            import traceback
            QMessageBox.critical(self,"PDF Hatasi",
                str(ex)+"\n\n"+traceback.format_exc()[-600:])


def make_pdf_multi(path, all_series, meta, flow, T, limit_pct=20, rsd_lim=5.0, lang="TR"):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.units import cm, mm
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                    Paragraph, Spacer, HRFlowable,
                                    PageBreak, KeepTogether)
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import os, math, numpy as np

    # Font
    fn = "Helvetica"; fb = "Helvetica-Bold"
    for nm, ff in [("DejaVu","DejaVuSans.ttf"),("DejaVuB","DejaVuSans-Bold.ttf")]:
        fp = resource_path(ff)
        if os.path.exists(fp):
            try: pdfmetrics.registerFont(TTFont(nm, fp))
            except: pass
    try:
        from reportlab.pdfbase.pdfmetrics import getFont
        getFont("DejaVuB"); fn="DejaVu"; fb="DejaVuB"
    except: pass

    W, H = A4
    BW = W - 3*cm   # tam sayfa genislik: 17.7 cm

    # Renksiz stiller
    def ps(sz, bold=False, color=colors.black, align=TA_LEFT):
        return ParagraphStyle("", fontName=fb if bold else fn,
            fontSize=sz, textColor=color, alignment=align,
            spaceBefore=1, spaceAfter=1, leading=sz+2)

    sTitle = ps(12, True,  colors.black, TA_CENTER)
    sSub   = ps(8,  False, colors.black, TA_CENTER)
    sHdr   = ps(8,  True,  colors.black, TA_LEFT)
    sLbl   = ps(7,  True,  colors.black, TA_LEFT)
    sVal   = ps(7,  False, colors.black, TA_CENTER)
    sValB  = ps(7,  True,  colors.black, TA_CENTER)
    sRed   = ps(7,  True,  colors.red,   TA_CENTER)

    # Tamamen renksiz tablo stili - sadece cizgiler
    def ts_clean(bold_header=True):
        return TableStyle([
            # Hic dolgu yok
            ("BACKGROUND", (0,0), (-1,-1), colors.white),
            ("FONTSIZE",   (0,0), (-1,-1), 7),
            ("FONTNAME",   (0,0), (-1,0),  fb),  # baslik bold
            ("ALIGN",      (0,0), (-1,-1), "CENTER"),
            ("ALIGN",      (0,0), (0,-1),  "LEFT"),   # ilk kolon sol
            ("VALIGN",     (0,0), (-1,-1), "MIDDLE"),
            ("TOPPADDING", (0,0), (-1,-1), 2),
            ("BOTTOMPADDING",(0,0),(-1,-1), 2),
            ("LEFTPADDING", (0,0),(-1,-1), 3),
            ("RIGHTPADDING",(0,0),(-1,-1), 3),
            # Dis cerceve kalin
            ("BOX",        (0,0), (-1,-1), 1.0, colors.black),
            # Baslik alti kalin cizgi
            ("LINEBELOW",  (0,0), (-1,0),  1.0, colors.black),
            # Ic yatay cizgiler ince
            ("INNERGRID",  (0,0), (-1,-1), 0.3, colors.black),
            # Alternatif satir golgesi yok - tamamen beyaz
        ])

    def make_section_title(text):
        """Bolum basligi - dolgu yok, sadece bold + alt cizgi"""
        return [
            Paragraph(text, ps(9, True, colors.black)),
            HRFlowable(width="100%", thickness=1.0, color=colors.black,
                       spaceAfter=2),
        ]

    def fmt(v, decimals=4):
        """Sayi formatlama: nokta yerine virgul"""
        if v is None: return "-"
        if isinstance(v, int): return str(v)
        try: return f"{v:.{decimals}f}".replace(".", ",")
        except: return str(v)

    doc = SimpleDocTemplate(path, pagesize=A4,
        leftMargin=1.5*cm, rightMargin=1.5*cm,
        topMargin=1.5*cm,  bottomMargin=1.5*cm,
        title="NGI Analysis Report")

    co = NGI_CUTOFFS[flow]
    story = []

    # =========================================================================
    # BASLIK BLOGU
    # =========================================================================
    if lang=="TR":
        pdf_main_title = "NGI Kaskadit Impaktor Analiz Araci"
        pdf_sub_title  = "Ph.Eur 2.9.18 / USP &lt;601&gt;  |  NGI Kaskadi Impaktoru Analizi"
    else:
        pdf_main_title = "Results and Analysis for Next Generation Impactor"
        pdf_sub_title  = "Ph.Eur 2.9.18 / USP &lt;601&gt;  |  NGI Cascade Impactor Analysis"
    story.append(Paragraph(pdf_main_title, sTitle))
    story.append(Paragraph(pdf_sub_title, sSub))
    story.append(HRFlowable(width="100%", thickness=1.5, color=colors.black,
                             spaceBefore=3, spaceAfter=3))

    # Meta bilgiler - tam genislik 2 kolonlu tablo
    meta_data = [
        [Paragraph("Date of Analysis:", sLbl),
         Paragraph(meta.get("date",""), sVal),
         Paragraph("Analyst:", sLbl),
         Paragraph(meta.get("operator",""), sVal),
         Paragraph("Flow Rate:", sLbl),
         Paragraph(f"{flow} L/min", sValB)],
        [Paragraph("Product Name:", sLbl),
         Paragraph(meta.get("product",""), sVal),
         Paragraph("Batch Number:", sLbl),
         Paragraph(meta.get("batch",""), sVal),
         Paragraph("Method:", sLbl),
         Paragraph("EP / Ph.Eur 2.9.18", sVal)],
    ]
    cw_meta = [2.5*cm, (BW-5*cm)/2-2.5*cm, 2.5*cm, (BW-5*cm)/2-2.5*cm, 2.5*cm, 2.5*cm]
    # 3 cift: label+value, esit pay
    cw_meta = [2.4*cm, 3.4*cm, 2.4*cm, 3.4*cm, 2.2*cm, 3.3*cm]
    mt = Table(meta_data, colWidths=cw_meta)
    mt.setStyle(TableStyle([
        ("BACKGROUND", (0,0),(-1,-1), colors.white),
        ("FONTSIZE",   (0,0),(-1,-1), 7),
        ("ALIGN",      (0,0),(-1,-1), "LEFT"),
        ("VALIGN",     (0,0),(-1,-1), "MIDDLE"),
        ("TOPPADDING", (0,0),(-1,-1), 2),
        ("BOTTOMPADDING",(0,0),(-1,-1), 2),
        ("LEFTPADDING", (0,0),(-1,-1), 3),
        ("BOX",        (0,0),(-1,-1), 0.8, colors.black),
        ("INNERGRID",  (0,0),(-1,-1), 0.3, colors.black),
    ]))
    story.append(mt)
    story.append(Spacer(1, 3*mm))

    # Cut-off tablosu - tam sayfa genisliginde
    vis_stages = [s for s in ["S1","S2","S3","S4","S5","S6","S7","MOC"]
                  if co.get(s,999) < 900]
    n_co = len(vis_stages) + 1
    co_cw = [BW/n_co] * n_co
    co_hdr = [Paragraph("Cut-off D50 (um)", sLbl)] + \
             [Paragraph(s, sLbl) for s in vis_stages]
    co_val = [Paragraph(f"{flow} L/min", sValB)] + \
             [Paragraph(fmt(co[s],3), sVal) for s in vis_stages]
    co_tbl = Table([co_hdr, co_val], colWidths=co_cw)
    co_tbl.setStyle(ts_clean())
    story.append(co_tbl)
    story.append(Spacer(1, 4*mm))

    # =========================================================================
    # HER SERI - HER RUN
    # =========================================================================
    for sd in all_series:
        ref_tag = "  [REFERENCE]" if sd["is_ref"] else ""
        for item in make_section_title(f"Series: {sd['name']}{ref_tag}"):
            story.append(item)
        story.append(Spacer(1, 2*mm))

        for run in sd["runs"]:
            run_items = []
            run_items.append(Paragraph(
                f"Run Number = {run['run_no']}  |  Sampling Flow Rate = {flow} L/min",
                ps(8, True, colors.black)))
            run_items.append(HRFlowable(width="100%", thickness=0.5,
                color=colors.black, spaceBefore=1, spaceAfter=2))

            if "error" in run:
                run_items.append(Paragraph(
                    f"Insufficient data (n={run.get('n',0)})", sVal))
                story.append(KeepTogether(run_items))
                story.append(Spacer(1,3*mm)); continue

            # Kumulatif tablo - tam sayfa genisliginde
            cum_cols = (["Kesme D50\n[um]","Kutle\n[mg]","Kumulatif\n[mg]",
                        "Kumulatif\n[%]","Gecerli","Probit z"]
                if lang=="TR" else
                ["Cut-off D50\n[um]","Mass\n[mg]","Cumulative mass\n[mg]",
                 "Cumulative\n[%]","Valid","Probit z"])
            n_cum = len(cum_cols)
            # Genislikler: D50=2cm, Mass=2.5cm, CumMass=2.5cm, Cum%=2.5cm, Valid=1.3cm, Probit=kalan
            cw_cum_fixed = [2.0*cm, 2.5*cm, 2.5*cm, 2.3*cm, 1.3*cm]
            cw_cum_last  = BW - sum(cw_cum_fixed)
            cw_cum = cw_cum_fixed + [cw_cum_last]

            cum_hdr_row = [Paragraph(h, sLbl) for h in cum_cols]
            cum_data = [cum_hdr_row]
            valid_st = {v["stage"] for v in run["valid"]}
            cum_m = 0.0
            ts_cum = ts_clean()
            row_idx = 1

            for row in run["cum_data"]:
                s = row["stage"]
                # Sadece ISM stage ve Throat goster
                if s not in (["Throat"] + [x for x in ALL_KEYS if co.get(x,999)<900]):
                    continue
                cum_m += row["mass"]
                iv = s in valid_st
                pz = ""
                if 0 < row["u_pct"] < 100:
                    try: pz = fmt(norm.ppf(row["u_pct"]/100), 4)
                    except: pass
                d50_str = fmt(row["d50"],3) if row["d50"]<900 else "---"
                cum_data.append([
                    Paragraph(d50_str, sVal),
                    Paragraph(fmt(row["mass"],4), sValB if iv else sVal),
                    Paragraph(fmt(cum_m,4), sVal),
                    Paragraph(fmt(row["u_pct"],3), sValB if iv else sVal),
                    Paragraph("*" if iv else "", sValB),
                    Paragraph(pz, sVal),
                ])
                if iv:
                    # Valid satirlari: sadece bold, dolgu yok
                    ts_cum.add("FONTNAME",(0,row_idx),(-1,row_idx), fb)
                    ts_cum.add("LINEABOVE",(0,row_idx),(-1,row_idx), 0.5, colors.black)
                    ts_cum.add("LINEBELOW",(0,row_idx),(-1,row_idx), 0.5, colors.black)
                row_idx += 1

            t_cum = Table(cum_data, colWidths=cw_cum, repeatRows=1)
            t_cum.setStyle(ts_cum)
            run_items.append(t_cum)
            run_items.append(Spacer(1,2*mm))

            # Parametreler - 2 satir, tam sayfa
            # Satir 1: FPD, MMAD, GSD, Intercept, Slope, R2, n
            p1_h = ["FPD (<=5um)\n[mg]","MMAD\n[um]","GSD",
                    "Intercept","Slope","R^2","n"]
            p1_v = [fmt(run["fpd"],4), fmt(run["mmad"],4), fmt(run["gsd"],4),
                    fmt(run["intercept"],3), fmt(run["slope"],3),
                    fmt(run["r2"],4), str(run["n"])]
            # Satir 2: FPF, Metered, Delivered
            p2_h = ["Fine Particle\nFraction [%]","Metered\n[mg]",
                    "Delivered\n[mg]","","","",""]
            p2_v = [fmt(run["fpf"],3), fmt(run["metered"],4),
                    fmt(run["delivered"],4),"","","",""]

            n_p = len(p1_h)
            cw_p = [BW/n_p]*n_p
            t_p = Table(
                [[Paragraph(h,sLbl) for h in p1_h],
                 [Paragraph(v,sValB if i<3 else sVal) for i,v in enumerate(p1_v)],
                 [Paragraph(h,sLbl) for h in p2_h],
                 [Paragraph(v,sValB if i<3 else sVal) for i,v in enumerate(p2_v)]],
                colWidths=cw_p)
            ts_p = ts_clean()
            # Baslik satirlari (0 ve 2) bold, deger satirlari (1 ve 3) normal
            ts_p.add("FONTNAME",(0,2),(-1,2), fb)
            ts_p.add("LINEABOVE",(0,2),(-1,2), 1.0, colors.black)
            t_p.setStyle(ts_p)
            run_items.append(t_p)
            run_items.append(Spacer(1,3*mm))
            story.append(KeepTogether(run_items))

        story.append(Spacer(1,2*mm))

    # =========================================================================
    # TABULAR SUMMARY
    # =========================================================================
    story.append(PageBreak())
    for item in make_section_title("Tabular Summary — NGI Stage Masses per Discharge"):
        story.append(item)
    story.append(Spacer(1,3*mm))

    for sd in all_series:
        valid_runs = [r for r in sd["runs"] if "error" not in r]
        if not valid_runs: continue
        n_r = len(valid_runs)

        story.append(Paragraph(
            f"Series: {sd['name']}" + ("  [REF]" if sd["is_ref"] else ""),
            ps(8, True, colors.black)))
        story.append(Spacer(1,1*mm))

        # Stage kutleleri - tam genislik
        disp_s = [s for s in ["Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]]
        ts_hdr = ([Paragraph("Stage",sLbl)] +
                  [Paragraph(f"Run {r['run_no']}",sLbl) for r in valid_runs] +
                  [Paragraph("Mean",sLbl),Paragraph("SD",sLbl),Paragraph("%RSD",sLbl)])
        n_ts_col = 1 + n_r + 3
        cw_ts_lbl = 1.8*cm
        cw_ts_rest = (BW - cw_ts_lbl) / (n_r+3)
        cw_ts = [cw_ts_lbl] + [cw_ts_rest]*(n_r+3)

        ts_rows = [ts_hdr]
        for s in disp_s:
            vals = [r["masses"].get(s,0) for r in valid_runs]
            if all(v==0 for v in vals): continue
            mean_v = float(np.mean(vals))
            sd_v   = float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            rsd_v  = sd_v/mean_v*100 if mean_v else 0.0
            ts_rows.append(
                [Paragraph(s, sLbl)] +
                [Paragraph(fmt(v,4),sVal) for v in vals] +
                [Paragraph(fmt(mean_v,4),sValB),
                 Paragraph(fmt(sd_v,4),sVal),
                 Paragraph(fmt(rsd_v,2),sVal)])
        t_ts = Table(ts_rows, colWidths=cw_ts, repeatRows=1)
        t_ts.setStyle(ts_clean())
        story.append(t_ts)
        story.append(Spacer(1,3*mm))

        # Parametre ozet - tam genislik
        story.append(Paragraph("Dose and Particle Size Characterisation",
            ps(7,True,colors.black)))
        story.append(Spacer(1,1*mm))
        pkeys = [("metered","Metered [mg]"),("delivered","Delivered [mg]"),
                 ("fpd","FPD 5µm [mg]"),("fpf","FPF 5µm [%]"),
                 ("fpd3","FPD 3µm [mg]"),("fpf3","FPF 3µm [%]"),
                 ("fpd15","FPD 1.5µm [mg]"),("fpf15","FPF 1.5µm [%]"),
                 ("mmad","MMAD [um]"),("gsd","GSD"),
                 ("d10","D10 [um]"),("d50","D50 [um]"),("d90","D90 [um]"),
                 ("slope","Slope"),("intercept","Intercept"),("r2","R^2")]
        p_hdr2 = ([Paragraph("Parameter",sLbl)] +
                  [Paragraph(f"Run {r['run_no']}",sLbl) for r in valid_runs] +
                  [Paragraph("Mean",sLbl),Paragraph("SD",sLbl),Paragraph("%RSD",sLbl)])
        cw_p2 = [2.8*cm] + [(BW-2.8*cm)/(n_r+3)]*(n_r+3)
        p_rows2 = [p_hdr2]
        for k,lbl in pkeys:
            vals2 = [r.get(k,0) for r in valid_runs if k in r]
            if not vals2: continue
            mean2=float(np.mean(vals2))
            sd2=float(np.std(vals2,ddof=1)) if len(vals2)>1 else 0.0
            rsd2=sd2/mean2*100 if mean2 else 0.0
            is_key = k in ("fpd","fpf","fpd3","fpf3","fpd15","fpf15","mmad","gsd","d10","d50","d90")
            sv = ps(7,True,colors.black,TA_CENTER) if is_key else sVal
            rsd_style = sRed if rsd2 > rsd_lim else sVal
            p_rows2.append(
                [Paragraph(lbl, sLbl if not is_key else ps(7,True,colors.black))] +
                [Paragraph(fmt(r.get(k,0),4),sv) for r in valid_runs] +
                [Paragraph(fmt(mean2,4),sv),
                 Paragraph(fmt(sd2,4),sVal),
                 Paragraph(fmt(rsd2,2),rsd_style)])
        t_p2 = Table(p_rows2, colWidths=cw_p2, repeatRows=1)
        t_p2.setStyle(ts_clean())
        story.append(t_p2)
        story.append(Spacer(1,5*mm))

    # =========================================================================
    # REFERANS KARSILASTIRMA
    # =========================================================================
    ref_sd = next((sd for sd in all_series if sd["is_ref"]), None)
    if ref_sd and ref_sd["avg"]:
        story.append(PageBreak())
        for item in make_section_title(
            f"Reference Comparison — {ref_sd['name']}  |  Limit: +/-{limit_pct:.0f}%"):
            story.append(item)
        story.append(Spacer(1,3*mm))

        ref_m = ref_sd["avg"]["avg_masses"]
        vis_cmp = ["Throat"] + [s for s in GRAPH_STAGES if co.get(s,999)<900]
        test_series = [sd for sd in all_series if not sd["is_ref"] and sd["avg"]]

        for sd in test_series:
            f2 = calc_f2(ref_m, sd["avg"]["avg_masses"], co)
            f2_pass = f2 is not None and f2 >= 50
            if f2 is not None:
                f2_txt = "f2 = " + fmt(f2,1) + "  " + \
                         ("Similar (>=50)" if f2_pass else "Different (<50)")
            else:
                f2_txt = "f2 = N/A"
            story.append(Paragraph(f"{sd['name']}   |   {f2_txt}",
                ps(8, True, colors.black)))
            story.append(Spacer(1,1*mm))

            cmp_hdr = [Paragraph(h,sLbl) for h in
                ["Stage","Reference\n[mg]","Test\n[mg]",
                 "Diff %",f"Upper\n+{limit_pct:.0f}%",
                 f"Lower\n-{limit_pct:.0f}%","Result"]]
            cmp_cw = [BW/7]*7
            cmp_rows = [cmp_hdr]
            ts_cmp = ts_clean()
            for ri2, s in enumerate(vis_cmp):
                rv = ref_m.get(s,0); tv = sd["avg"]["avg_masses"].get(s,0)
                diff = (tv-rv)/rv*100 if rv else 0.0
                in_lim = abs(diff) <= limit_pct
                res_clr = colors.HexColor("#006600") if in_lim else colors.red
                cmp_rows.append([
                    Paragraph(s, sLbl),
                    Paragraph(fmt(rv,4), sVal),
                    Paragraph(fmt(tv,4), sVal),
                    Paragraph(f"{diff:+,.2f}%".replace('.',','), sVal),
                    Paragraph(fmt(rv*(1+limit_pct/100),4), sVal),
                    Paragraph(fmt(rv*(1-limit_pct/100),4), sVal),
                    Paragraph("PASS" if in_lim else "FAIL",
                        ps(7,True,colors.black,TA_CENTER)),
                ])
                if not in_lim:
                    ts_cmp.add("FONTNAME",(0,ri2+1),(-1,ri2+1),fb)
            t_cmp = Table(cmp_rows, colWidths=cmp_cw, repeatRows=1)
            t_cmp.setStyle(ts_cmp)
            story.append(t_cmp)
            story.append(Spacer(1,4*mm))

    # =========================================================================
    # GRAFIK SAYFASI
    # =========================================================================
    story.append(PageBreak())
    for item in make_section_title("Graphical Analysis"):
        story.append(item)
    story.append(Spacer(1,3*mm))

    import io as _io
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    from reportlab.platypus import Image as RLImage

    BW_COLORS = ["black","#333","#555","#777","#999","#BBB"]
    BW_LS     = ["-","--","-.",":","--","-."]

    # Log-Probit
    fig_lp, ax_lp = _plt.subplots(figsize=(7,3.8))
    ax_lp.set_facecolor("white"); fig_lp.patch.set_facecolor("white")
    for si, sd in enumerate(all_series):
        valid_runs = [r for r in sd["runs"] if "error" not in r]
        if not valid_runs: continue
        col = BW_COLORS[si % len(BW_COLORS)]
        ls  = BW_LS[si % len(BW_LS)]
        lw  = 2.0 if sd["is_ref"] else 1.3
        min_len = min(len(r["x_reg"]) for r in valid_runs)
        avg_x = sum(r["x_reg"][:min_len] for r in valid_runs)/len(valid_runs)
        avg_y = sum(r["y_reg"][:min_len] for r in valid_runs)/len(valid_runs)
        b_avg = sum(r["b"] for r in valid_runs)/len(valid_runs)
        a_avg = sum(r["a"] for r in valid_runs)/len(valid_runs)
        ax_lp.plot(avg_x, avg_y, "o", color=col, ms=5, zorder=4)
        xr = np.linspace(min(avg_x)-0.15, max(avg_x)+0.15, 60)
        ax_lp.plot(xr, a_avg+b_avg*xr, ls=ls, color=col, lw=lw,
            label=sd["name"]+(" [REF]" if sd["is_ref"] else ""))
        if sd["avg"] and "mmad" in sd["avg"]["params"]:
            mv = sd["avg"]["params"]["mmad"][0]
            if mv > 0:
                ax_lp.axvline(math.log10(mv), color=col, lw=0.7, ls=":", alpha=0.8)
    ax_lp.set_xlabel("log10(D50, um)", fontsize=9)
    ax_lp.set_ylabel("Probit z", fontsize=9)
    ax_lp.set_title(f"Log-Probit  [{flow} L/min]  (Series Averages)", fontsize=10)
    ax_lp.legend(fontsize=7, framealpha=0.9)
    ax_lp.grid(True, ls="--", alpha=0.3, color="gray")
    ax_lp.tick_params(labelsize=8)
    for sp in ax_lp.spines.values(): sp.set_color("black"); sp.set_linewidth(0.8)
    fig_lp.tight_layout()
    buf_lp = _io.BytesIO()
    fig_lp.savefig(buf_lp, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf_lp.seek(0); _plt.close(fig_lp)
    story.append(RLImage(buf_lp, width=BW, height=BW*3.8/7))

    # Log-Probit alti: MMAD degerleri tablosu
    mmad_rows = [[Paragraph("Series", sLbl),
                  Paragraph("MMAD [um]", sLbl),
                  Paragraph("D10 [um]", sLbl),
                  Paragraph("D50 [um]", sLbl),
                  Paragraph("D90 [um]", sLbl),
                  Paragraph("GSD", sLbl),
                  Paragraph("Slope", sLbl),
                  Paragraph("Intercept", sLbl),
                  Paragraph("R^2", sLbl)]]
    for sd2 in all_series:
        if not sd2["avg"]: continue
        pr = sd2["avg"]["params"]
        mmad_rows.append([
            Paragraph(sd2["name"] + (" [REF]" if sd2["is_ref"] else ""), sLbl),
            Paragraph(fmt(pr.get("mmad",(0,))[0], 3), sValB),
            Paragraph(fmt(pr.get("gsd",(0,))[0], 3), sVal),
            Paragraph(fmt(pr.get("slope",(0,))[0], 3), sVal),
            Paragraph(fmt(pr.get("intercept",(0,))[0], 3), sVal),
            Paragraph(fmt(pr.get("r2",(0,))[0], 4), sVal),
        ])
    n_mc = len(mmad_rows[0])
    t_mmad = Table(mmad_rows, colWidths=[BW/n_mc]*n_mc)
    t_mmad.setStyle(ts_clean())
    story.append(t_mmad)
    story.append(Spacer(1,5*mm))

    # APSD
    vis_ap = ["Throat"] + [s for s in ["S1","S2","S3","S4","S5","S6","S7","MOC"]
                            if co.get(s,999)<900]
    x_ap = list(range(len(vis_ap)))
    fig_ap, ax_ap = _plt.subplots(figsize=(7,3.8))
    ax_ap.set_facecolor("white"); fig_ap.patch.set_facecolor("white")
    ref_ap = None
    for si, sd in enumerate(all_series):
        if not sd["avg"]: continue
        col = BW_COLORS[si % len(BW_COLORS)]
        ls  = BW_LS[si % len(BW_LS)]
        lw  = 2.0 if sd["is_ref"] else 1.3
        ms_ap = [sd["avg"]["avg_masses"].get(s,0) for s in vis_ap]
        valid_runs = [r for r in sd["runs"] if "error" not in r]
        sds_ap = []
        for s in vis_ap:
            vals=[r["masses"].get(s,0) for r in valid_runs]
            sds_ap.append(float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0)
        ax_ap.plot(x_ap, ms_ap, color=col, lw=lw, ms=5, marker="o", linestyle=ls if isinstance(ls,str) else "-",
            label=sd["name"]+(" [REF]" if sd["is_ref"] else ""))
        ax_ap.fill_between(x_ap, [m-s for m,s in zip(ms_ap,sds_ap)],
            [m+s for m,s in zip(ms_ap,sds_ap)], color=col, alpha=0.08)
        if sd["is_ref"]: ref_ap = ms_ap
    if ref_ap:
        up = [v*(1+limit_pct/100) for v in ref_ap]
        lo = [v*(1-limit_pct/100) for v in ref_ap]
        ax_ap.plot(x_ap, up, "--", color="black", lw=0.8, alpha=0.7,
            label=f"+{limit_pct:.0f}%")
        ax_ap.plot(x_ap, lo, "--", color="black", lw=0.8, alpha=0.7,
            label=f"-{limit_pct:.0f}%")
    ax_ap.set_xticks(x_ap); ax_ap.set_xticklabels(vis_ap, rotation=20, ha="right", fontsize=8)
    ax_ap.set_xlabel("Stage", fontsize=9)
    ax_ap.set_ylabel("Mean Mass (mg/actuation)", fontsize=9)
    ax_ap.set_title(f"APSD Distribution  [{flow} L/min]  (Series Averages +/- SD)", fontsize=10)
    ax_ap.legend(fontsize=7, framealpha=0.9)
    ax_ap.grid(True, ls="--", alpha=0.3, color="gray")
    ax_ap.tick_params(labelsize=8)
    for sp in ax_ap.spines.values(): sp.set_color("black"); sp.set_linewidth(0.8)
    fig_ap.tight_layout()
    buf_ap = _io.BytesIO()
    fig_ap.savefig(buf_ap, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf_ap.seek(0); _plt.close(fig_ap)
    story.append(RLImage(buf_ap, width=BW, height=BW*3.8/7))

    # Footer
    story.append(Spacer(1,4*mm))
    story.append(HRFlowable(width="100%", thickness=0.8, color=colors.black))
    story.append(Paragraph(
        f"NGI Analysis Tool v5  |  Ph.Eur 2.9.18 / USP &lt;601&gt;  |  "
        f"Generated: {datetime.now().strftime('%d.%m.%Y %H:%M')}",
        ps(6.5, False, colors.black, TA_CENTER)))

    doc.build(story)



if __name__ == "__main__":
    app=QApplication(sys.argv)
    app.setApplicationName("NGI Impactor Analysis")
    win=NGIApp()
    win.show()
    sys.exit(app.exec())
