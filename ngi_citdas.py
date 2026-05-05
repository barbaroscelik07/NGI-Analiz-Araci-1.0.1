"""NGI Cascade Impactor Analysis Tool v5 - CITDAS validated"""
import tkinter as tk
from tkinter import filedialog, messagebox
import customtkinter as ctk
import numpy as np
from scipy.stats import norm
import math, os, sys
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from datetime import datetime

def resource_path(rel):
    base = getattr(sys,'_MEIPASS', os.path.dirname(os.path.abspath(
        sys.executable if getattr(sys,'frozen',False) else __file__)))
    return os.path.join(base, rel)

NGI_CUTOFFS = {
    15: {"Device":999,"Throat":999,"Presep":999,"S1":999,
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
STAGE_ORDER   = ["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7"]
ALL_KEYS      = STAGE_ORDER + ["MOC"]
ISM_STAGES    = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
GRAPH_STAGES  = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
RUNS_PER_SERIES = 3
EXCLUSIVE_FLOWS = {15}
CP = ["#2E75B6","#ED7D31","#70AD47","#E84040","#7030A0",
      "#00B0F0","#D4A000","#C00000","#00B050","#FF69B4"]

L = {
"TR":{
 "title":"NGI Impaktor Analiz Araci",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | Next Generation Impactor",
 "lang_btn":"English","product":"Urun Adi","batch":"Lot No.",
 "operator":"Analist","date":"Tarih","flow_rate":"Akis Hizi",
 "add_series":"+ Seri Ekle","del_series":"Seri Sil",
 "calculate":"Hesapla","clear":"Temizle","export_pdf":"PDF Rapor",
 "tab_results":"Sonuclar","tab_plot":"Log-Probit",
 "tab_dist":"Dagilim","tab_summary":"Ozet","tab_compare":"Karsilastirma",
 "series":"Seri","run":"Run","paste_btn":"Yapistir",
 "mean":"Ort.","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Yetersiz nokta","status_ready":"Hazir.",
 "status_done":"Hesaplama tamamlandi.",
 "valid_range":"Gecerlilik (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"Bu seri REFERANS","ref_label":"REFERANS",
 "limit_label":"Limit Tipi","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manuel (%)","lim_pct":"Limit %",
 "f2_label":"f2 Benzerlik Faktoru","f2_pass":">=50 Benzer","f2_fail":"<50 Farkli",
 "outside_warn":"UYARI: Limit disi noktalar","no_ref":"Referans secilmedi",
 "ddu_label":"DDU Analizi","trend_label":"Trend Grafigi",
 "rsd_limit":"RSD Kabul (%)","cv_label":"CV%","load_csv":"CSV Yukle","csv_template":"Sablon","csv_loaded":"CSV yuklendi: {n} seri, {r} run","csv_err_format":"CSV format hatasi: Beklenen kolonlar eksik","csv_err_flow":"Uyari: Farkli flow hizlari tespit edildi","csv_ref_ask":"Referans seri secin","csv_ref_none":"Referans yok (atla)","csv_4runs":"Uyari: {s} serisi 3den fazla run iceriyor, ilk 3 alindi",
 "stage":"Stage","mass_mg":"Kutle (mg)","cum_mass":"Kum. Kutle",
 "cum_pct":"Kum. %","valid_pt":"Gecerli","probit_z":"Probit z",
 "param":"Parametre","accept":"Kabul","fail_lbl":"BASARISIZ",
 "pass_lbl":"GECTI","ddu_mean":"Ort.","ddu_sd":"SD","ddu_rsd":"RSD%",
 "fp_dose":"FPD (mg)","fp_frac":"FPF (%)","slope_lbl":"Egim (Slope)",
 "int_lbl":"Intercept","r2_lbl":"R2","n_lbl":"n",
 "ref_comp":"Referans Karsilastirma","series_avg":"Seri Ortalamasi",
 "trend_mmad":"MMAD Trendi","trend_gsd":"GSD Trendi",
 "dec_sep":","
},
"EN":{
 "title":"NGI Cascade Impactor Analysis",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | Next Generation Impactor",
 "lang_btn":"Turkce","product":"Product","batch":"Batch No.",
 "operator":"Analyst","date":"Date","flow_rate":"Flow Rate",
 "add_series":"+ Add Series","del_series":"Del Series",
 "calculate":"Calculate","clear":"Clear","export_pdf":"PDF Report",
 "tab_results":"Results","tab_plot":"Log-Probit",
 "tab_dist":"Distribution","tab_summary":"Summary","tab_compare":"Compare",
 "series":"Series","run":"Run","paste_btn":"Paste",
 "mean":"Mean","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Insufficient pts","status_ready":"Ready.",
 "status_done":"Calculation complete.",
 "valid_range":"Valid Range (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"This series is REFERENCE","ref_label":"REFERENCE",
 "limit_label":"Limit Type","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manual (%)","lim_pct":"Limit %",
 "f2_label":"f2 Similarity Factor","f2_pass":">=50 Similar","f2_fail":"<50 Different",
 "outside_warn":"WARNING: Points outside limits","no_ref":"No reference selected",
 "ddu_label":"DDU Analysis","trend_label":"Trend Chart",
 "rsd_limit":"RSD Accept (%)","cv_label":"CV%","load_csv":"Load CSV","csv_template":"Template","csv_loaded":"CSV loaded: {n} series, {r} runs","csv_err_format":"CSV format error: Expected columns missing","csv_err_flow":"Warning: Multiple flow rates detected","csv_ref_ask":"Select reference series","csv_ref_none":"No reference (skip)","csv_4runs":"Warning: {s} has more than 3 runs, first 3 taken",
 "stage":"Stage","mass_mg":"Mass (mg)","cum_mass":"Cum. Mass",
 "cum_pct":"Cum. %","valid_pt":"Valid","probit_z":"Probit z",
 "param":"Parameter","accept":"Accept","fail_lbl":"FAIL",
 "pass_lbl":"PASS","ddu_mean":"Mean","ddu_sd":"SD","ddu_rsd":"RSD%",
 "fp_dose":"FPD (mg)","fp_frac":"FPF (%)","slope_lbl":"Slope",
 "int_lbl":"Intercept","r2_lbl":"R2","n_lbl":"n",
 "ref_comp":"Reference Comparison","series_avg":"Series Summary",
 "trend_mmad":"MMAD Trend","trend_gsd":"GSD Trend",
 "dec_sep":"."
}}

def calc_run(masses, flow, lo=15, hi=85, delivered_tp=False):
    """delivered_tp=True: Delivered = Throat+Presep+ISM, False: ISM only"""
    co   = NGI_CUTOFFS[flow]
    excl = flow in EXCLUSIVE_FLOWS
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
                any_v=any(lo<u<hi for _,u in [pts_all[i],pts_all[i+1]])
                if any_v:
                    t=(tu-u1)/(u2-u1)
                    return 10**(math.log10(d1)+t*(math.log10(d2)-math.log10(d1)))
        return 10**((tz-a)/b)
    d84=get_d(84.13,1.0); d16=get_d(15.87,-1.0)
    gsd=math.sqrt(d84/d16) if (d84 and d16 and d16>0) else 10**(1/b)
    d5u=None
    for i in range(len(pts_all)-1):
        d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
        if d1<=5<=d2 and d2>d1:
            t=(math.log10(5)-math.log10(d1))/(math.log10(d2)-math.log10(d1))
            d5u=u1+t*(u2-u1); break
    if d5u is None:
        d5u=norm.cdf(a+b*math.log10(5))*100
    fpd=d5u/100*ism; fpf=fpd/metered*100 if metered>0 else 0
    res.update({"n":len(valid),"a":a,"b":b,"slope":b,"intercept":a+5,"r2":r2,
                "mmad":mmad,"gsd":gsd,"fpd":fpd,"fpf":fpf,"x_reg":x,"y_reg":y})
    return res

def calc_series_avg(runs):
    valid=[r for r in runs if "error" not in r]
    if not valid: return None
    avg_masses={}
    for s in ALL_KEYS:
        vals=[r["masses"].get(s,0) for r in valid]
        avg_masses[s]=float(np.mean(vals))
    params={}
    for p in ["mmad","gsd","fpd","fpf","metered","delivered","slope","intercept","r2"]:
        vals=[r[p] for r in valid if p in r]
        if vals:
            m=float(np.mean(vals)); s=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            params[p]=(m,s,s/m*100 if m else 0.0)
    return {"avg_masses":avg_masses,"params":params,"n_valid":len(valid)}

def calc_f2(ref_m, test_m, co):
    stages=[s for s in GRAPH_STAGES if co.get(s,999)<900]
    diffs=[]
    for s in stages:
        r=ref_m.get(s,0); t=test_m.get(s,0)
        if r>0: diffs.append(((t-r)/r*100)**2)
    if not diffs: return None
    return 50*math.log10(100/math.sqrt(1+np.mean(diffs)))

def parse_paste(text):
    import re
    lines=[l.strip() for l in text.strip().splitlines() if l.strip()]
    if not lines: return None
    def split_line(line):
        # Sekme VEYA çoklu boşluk ayırıcı
        if '\t' in line:
            return [t.strip() for t in line.split('\t')]
        else:
            return [t.strip() for t in re.split(r'\s{2,}|\s+', line.strip())]
    def is_header(line):
        parts = split_line(line)
        first = parts[0] if parts else ''
        try: float(first.replace(',','.').replace(' ','')); return False
        except: return True
    if is_header(lines[0]): lines=lines[1:]
    if not lines: return None
    result=[]
    for line in lines:
        tokens=split_line(line)
        try:
            vals=[float(t.replace(',','.').replace(' ','')) for t in tokens if t]
            if len(vals)>=11: result.append(vals[:11])
            elif len(vals)>=10:  # Device eksik olabilir
                result.append([0.0]+vals[:10])
        except: pass
    return result if result else None

DISP_STAGES=["Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]

def parse_csv(path):
    """CSV dosyasini oku, seri/run gruplarini dondur.
    Return: (series_dict, flow, warnings)
    series_dict: {seri_adi: {"runs":[{stage:val}], "ref":bool}}
    """
    import csv, re
    warnings = []
    series_dict = {}  # {ad: {"runs":[], "ref":False}}
    flow_vals = set()

    with open(path, newline="", encoding="utf-8-sig") as f:
        raw = f.read()
    lines = [l.strip() for l in raw.splitlines() if l.strip()]
    if not lines:
        return None, None, ["Dosya bos"]

    # Ayrac tespit: tum dosyaya bak
    sep = ";" if raw.count(";") > raw.count(",") else ","

    # Gercek header satirini bul: "seri" VE ("run" veya "flow") iceren ilk satir
    # Basi aciklama satirlarini otomatik atla
    hdr_idx = None
    for i, line in enumerate(lines):
        cols = [c.strip().lower() for c in line.split(sep)]
        if "seri" in cols and ("run" in cols or "flow" in cols):
            hdr_idx = i
            break
    if hdr_idx is None:
        return None, None, ["Header satiri bulunamadi. 'Seri', 'Run', 'Flow' kolonlari olmali."]

    hdr = [h.strip().lower().replace(" ","_") for h in lines[hdr_idx].split(sep)]
    lines = lines[hdr_idx+1:]  # Veri satirlari header sonrasindan baslar

    # Zorunlu kolonlari kontrol et
    missing = [c for c in ["seri","run","flow","throat","s1"] if c not in hdr]
    if missing:
        return None, None, [f"Eksik kolon: {missing}"]

    def to_float(s):
        if s is None or str(s).strip() == "": return 0.0
        v = str(s).strip()
        # Hem 1.234 hem 1,234 formatini destekle
        # Eger sadece bir virgul varsa ve nokta yoksa: ondalik virgul
        if v.count(",") == 1 and "." not in v:
            v = v.replace(",", ".")
        elif v.count(",") > 1:
            # Binlik ayrac virgul: "1,234.56" -> "1234.56"
            v = v.replace(",", "")
        return float(v)

    ref_col = "referans" in hdr

    for line in lines[1:]:
        # Bos satiri atla
        if not line.replace(",","").strip():
            continue
        parts = line.split(sep)
        # Kolonlari map et
        row = {}
        for i, h in enumerate(hdr):
            row[h] = parts[i].strip() if i < len(parts) else ""

        seri = row.get("seri","").strip()
        if not seri: continue

        try:
            run_no = int(float(row.get("run","1") or "1"))
            flow_v = int(float(row.get("flow","60") or "60"))
        except:
            continue

        flow_vals.add(flow_v)

        is_ref = False
        if ref_col:
            ref_val = row.get("referans","0").strip()
            is_ref = ref_val in ("1","1.0","evet","yes","true","referans","ref")

        # Kutleleri al
        masses = {
            "Device": 0.0,
            "Throat": to_float(row.get("throat")),
            "Presep": to_float(row.get("presep")),
            "S1": to_float(row.get("s1")),
            "S2": to_float(row.get("s2")),
            "S3": to_float(row.get("s3")),
            "S4": to_float(row.get("s4")),
            "S5": to_float(row.get("s5")),
            "S6": to_float(row.get("s6")),
            "S7": to_float(row.get("s7")),
            "MOC": to_float(row.get("moc")),
        }

        if seri not in series_dict:
            series_dict[seri] = {"runs": [], "ref": is_ref, "flow": flow_v}
        if is_ref:
            series_dict[seri]["ref"] = True

        # Max 3 run
        if len(series_dict[seri]["runs"]) < 3:
            series_dict[seri]["runs"].append({"run_no": run_no, "masses": masses})
        else:
            warnings.append(f"csv_4runs__{seri}")

    if len(flow_vals) > 1:
        warnings.append("csv_err_flow")

    flow = list(flow_vals)[0] if flow_vals else 60
    return series_dict, flow, warnings


def fmt_num(v, decimals=4, dec_sep=","):
    """Sayiyi formatlayip ondalik ayracini uygula"""
    if isinstance(v, int): return str(v)
    s = f"{v:.{decimals}f}"
    return s.replace(".", dec_sep) if dec_sep != "." else s

class NGIApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.lang="TR"; self.T=L["TR"]
        self.all_series=[]; self.series_widgets=[]
        self.flow=60; self.lo=15; self.hi=85
        self.ref_var=tk.BooleanVar(value=False)
        self.var_lp_avg_only=tk.BooleanVar(value=False)
        self.limit_var=tk.StringVar(value="ema")
        self.custom_pct_var=tk.StringVar(value="20")
        self.rsd_limit_var=tk.StringVar(value="5")
        self.var_delivered_tp=tk.BooleanVar(value=False)
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")
        self.title(self.T["title"])
        self.geometry("1520x980"); self.minsize(1200,780)
        ico=resource_path("icon.ico")
        if os.path.exists(ico):
            try: self.iconbitmap(ico)
            except: pass
        self._build_ui()

    def _build_ui(self):
        hdr=ctk.CTkFrame(self,fg_color="#002D62",corner_radius=0,height=52)
        hdr.pack(fill="x"); hdr.pack_propagate(False)
        self.lbl_title=ctk.CTkLabel(hdr,text=self.T["title"],
            font=ctk.CTkFont(size=15,weight="bold"),text_color="#FFC600")
        self.lbl_title.pack(side="left",padx=14,pady=6)
        self.lbl_sub=ctk.CTkLabel(hdr,text=self.T["subtitle"],
            font=ctk.CTkFont(size=10),text_color="#aac8e8")
        self.lbl_sub.pack(side="left",padx=4)
        self.btn_lang=ctk.CTkButton(hdr,text=self.T["lang_btn"],width=80,height=28,
            command=self._toggle_lang,fg_color="#001a40",hover_color="#003580")
        self.btn_lang.pack(side="right",padx=12)
        body=ctk.CTkFrame(self,fg_color="transparent"); body.pack(fill="both",expand=True)
        self.left=ctk.CTkScrollableFrame(body,width=470,fg_color="#141824",corner_radius=0)
        self.left.pack(side="left",fill="y")
        self._build_left()
        right=ctk.CTkFrame(body,fg_color="#0e1219"); right.pack(side="left",fill="both",expand=True)
        self._build_right(right)
        sb=ctk.CTkFrame(self,height=24,fg_color="#090c12",corner_radius=0)
        sb.pack(fill="x",side="bottom"); sb.pack_propagate(False)
        self.lbl_status=ctk.CTkLabel(sb,text=self.T["status_ready"],
            anchor="w",font=ctk.CTkFont(size=10),text_color="#7090b0")
        self.lbl_status.pack(side="left",padx=8)

    def _build_left(self):
        p=self.left
        mf=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=8)
        mf.pack(fill="x",padx=6,pady=(8,4))
        for row_i,(key,default) in enumerate([
            ("product",""),("batch",""),("operator",""),
            ("date",datetime.now().strftime("%d.%m.%Y"))
        ]):
            ctk.CTkLabel(mf,text=self.T[key],font=ctk.CTkFont(size=11),
                width=82,anchor="e").grid(row=row_i,column=0,padx=(8,4),pady=3,sticky="e")
            e=ctk.CTkEntry(mf,height=28,font=ctk.CTkFont(size=11))
            e.insert(0,default); e.grid(row=row_i,column=1,padx=(0,8),pady=3,sticky="ew")
            setattr(self,f"e_{key}",e)
        mf.columnconfigure(1,weight=1)
        ff=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=8)
        ff.pack(fill="x",padx=6,pady=4)
        ctk.CTkLabel(ff,text=self.T["flow_rate"],
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFC600").pack(side="left",padx=(10,4),pady=6)
        self.var_flow=ctk.StringVar(value="60")
        ctk.CTkOptionMenu(ff,values=[str(f) for f in sorted(NGI_CUTOFFS.keys())],
            variable=self.var_flow,command=self._on_flow,width=70,height=28,
            font=ctk.CTkFont(size=12,weight="bold")).pack(side="left",padx=4)
        ctk.CTkLabel(ff,text="L/min",font=ctk.CTkFont(size=11),
            text_color="#aac8e8").pack(side="left",padx=(0,12))
        ctk.CTkLabel(ff,text=self.T["valid_range"],font=ctk.CTkFont(size=11)).pack(side="left")
        self.e_lo=ctk.CTkEntry(ff,width=48,height=28,justify="center",font=ctk.CTkFont(size=11))
        self.e_lo.insert(0,"15"); self.e_lo.pack(side="left",padx=3)
        ctk.CTkLabel(ff,text="-",font=ctk.CTkFont(size=11)).pack(side="left")
        self.e_hi=ctk.CTkEntry(ff,width=48,height=28,justify="center",font=ctk.CTkFont(size=11))
        self.e_hi.insert(0,"85"); self.e_hi.pack(side="left",padx=3)
        self.cbox=ctk.CTkFrame(p,fg_color="#111827",corner_radius=6)
        self.cbox.pack(fill="x",padx=6,pady=(2,4))
        self._refresh_cutoffs()
        # Birinci buton satiri: Seri + Hesap
        bf=ctk.CTkFrame(p,fg_color="transparent"); bf.pack(fill="x",padx=6,pady=(4,2))
        self.btn_add_s=ctk.CTkButton(bf,text=self.T["add_series"],width=110,height=30,
            command=self._add_series,fg_color="#1a3a6a",hover_color="#2255a0",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_add_s.pack(side="left",padx=(0,4))
        self.btn_del_s=ctk.CTkButton(bf,text=self.T["del_series"],width=90,height=30,
            command=self._del_series,fg_color="#3a1a1a",hover_color="#6a2020",
            font=ctk.CTkFont(size=11)); self.btn_del_s.pack(side="left",padx=(0,4))
        self.btn_calc=ctk.CTkButton(bf,text=self.T["calculate"],width=90,height=30,
            command=self._calculate,fg_color="#1a5a1a",hover_color="#2a8a2a",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_calc.pack(side="left",padx=(0,4))
        self.btn_clr=ctk.CTkButton(bf,text=self.T["clear"],width=70,height=30,
            command=self._clear,fg_color="#3a3a1a",hover_color="#6a6a20",
            font=ctk.CTkFont(size=11)); self.btn_clr.pack(side="left")
        # Ikinci buton satiri: PDF + CSV
        bf2=ctk.CTkFrame(p,fg_color="transparent"); bf2.pack(fill="x",padx=6,pady=(0,4))
        self.btn_pdf=ctk.CTkButton(bf2,text=self.T["export_pdf"],width=130,height=30,
            command=self._export_pdf,fg_color="#3a1a5a",hover_color="#6a20a0",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_pdf.pack(side="left",padx=(0,4))
        self.btn_csv=ctk.CTkButton(bf2,text=self.T["load_csv"],width=130,height=30,
            command=self._load_csv,fg_color="#0a4a3a",hover_color="#0a7a5a",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_csv.pack(side="left")
        rf2=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=6)
        rf2.pack(fill="x",padx=6,pady=(0,4))
        ctk.CTkLabel(rf2,text=self.T["rsd_limit"],font=ctk.CTkFont(size=11),
            text_color="#aac8e8").pack(side="left",padx=(8,4),pady=4)
        ctk.CTkEntry(rf2,textvariable=self.rsd_limit_var,width=48,height=24,
            justify="center",font=ctk.CTkFont(size=11)).pack(side="left",padx=4)
        ctk.CTkLabel(rf2,text="%  ",font=ctk.CTkFont(size=11)).pack(side="left")
        ctk.CTkCheckBox(rf2,text="Delivered = T+P+ISM",variable=self.var_delivered_tp,
            font=ctk.CTkFont(size=10),text_color="#aac8e8").pack(side="left",padx=(8,4))
        self.series_box=ctk.CTkFrame(p,fg_color="transparent")
        self.series_box.pack(fill="x",padx=4,pady=4)
        self._add_series()

    def _on_tab_change(self):
        """Sekme degisince sadece o sekmeyi render et"""
        if not self.all_series: return
        current = self.tabs.get()
        if current == self._last_rendered_tab: return
        self._last_rendered_tab = current
        T = self.T
        if current == T["tab_plot"]:    self._plot_lp()
        elif current == T["tab_dist"]:  self._plot_dist()
        elif current == T["tab_summary"]: self._show_summary()
        elif current == T["tab_compare"]: self._show_compare()
        # tab_results her zaman ilk render edilir (_calculate'de)

    def _refresh_cutoffs(self):
        for w in self.cbox.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        vis=[s for s in ["S2","S3","S4","S5","S6","S7","MOC"] if co.get(s,999)<900]
        hf=ctk.CTkFrame(self.cbox,fg_color="transparent"); hf.pack(fill="x",padx=4,pady=(3,0))
        ctk.CTkLabel(hf,text=f"{self.T['cutoff_title']}  [{flow} L/min]",
            font=ctk.CTkFont(size=10,weight="bold"),text_color="#FFC600").pack(side="left")
        vf=ctk.CTkFrame(self.cbox,fg_color="transparent"); vf.pack(fill="x",padx=4,pady=(1,4))
        for s in vis:
            sf=ctk.CTkFrame(vf,fg_color="#1a2540",corner_radius=4); sf.pack(side="left",padx=2)
            ctk.CTkLabel(sf,text=s,font=ctk.CTkFont(size=9,weight="bold"),
                text_color="#7ab0d0",width=30).pack(pady=(1,0))
            ctk.CTkLabel(sf,text=f"{co[s]:.2f}",font=ctk.CTkFont(size=9),
                text_color="#e0f0ff",width=30).pack(pady=(0,1))

    def _on_flow(self,v): self.flow=int(v); self._refresh_cutoffs()

    def _add_series(self):
        idx=len(self.series_widgets)+1; color=CP[(idx-1)%len(CP)]
        frame=ctk.CTkFrame(self.series_box,fg_color="#1c2336",corner_radius=8)
        frame.pack(fill="x",pady=4,padx=2)
        hf=ctk.CTkFrame(frame,fg_color="#001a40",corner_radius=6); hf.pack(fill="x",padx=6,pady=(6,2))
        ctk.CTkLabel(hf,text=f"  {self.T['series']} {idx}",
            font=ctk.CTkFont(size=12,weight="bold"),text_color=color).pack(side="left",pady=4)
        name_var=ctk.StringVar(value=f"Seri {idx}")
        ctk.CTkEntry(hf,textvariable=name_var,height=26,width=120,
            font=ctk.CTkFont(size=11)).pack(side="left",padx=8)
        paste_btn=ctk.CTkButton(hf,text=self.T["paste_btn"],width=80,height=26,
            font=ctk.CTkFont(size=10),fg_color="#003580",hover_color="#0055c0")
        paste_btn.pack(side="right",padx=6)
        ref_check=None
        if idx==1:
            rf=ctk.CTkFrame(frame,fg_color="transparent"); rf.pack(fill="x",padx=6,pady=(0,2))
            ref_check=ctk.CTkCheckBox(rf,text=self.T["ref_check"],variable=self.ref_var,
                font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFD700",
                fg_color="#8B6914",hover_color="#B8860B",border_color="#FFD700")
            ref_check.pack(side="left",padx=4)
        grid=ctk.CTkFrame(frame,fg_color="transparent"); grid.pack(fill="x",padx=8,pady=(2,6))
        ctk.CTkLabel(grid,text="Stage",width=64,font=ctk.CTkFont(size=11,weight="bold"),
            text_color="#5a8ab0",anchor="w").grid(row=0,column=0,padx=2,pady=1)
        for ri in range(RUNS_PER_SERIES):
            ctk.CTkLabel(grid,text=f"Run {ri+1}",width=78,
                font=ctk.CTkFont(size=11,weight="bold"),
                text_color=color,anchor="center").grid(row=0,column=ri+1,padx=2,pady=1)
        run_entries=[{} for _ in range(RUNS_PER_SERIES)]
        for si,s in enumerate(DISP_STAGES):
            row_i=si+1
            lc="#FFD700" if s=="Presep" else "#aac8e8"
            ctk.CTkLabel(grid,text=s,width=64,font=ctk.CTkFont(size=11),
                text_color=lc,anchor="w").grid(row=row_i,column=0,padx=2,pady=1)
            for ri in range(RUNS_PER_SERIES):
                v=ctk.StringVar(value="0.000")
                e=ctk.CTkEntry(grid,textvariable=v,height=24,width=78,
                    font=ctk.CTkFont(size=11),justify="center")
                e.grid(row=row_i,column=ri+1,padx=2,pady=1)
                e.bind("<FocusIn>",  lambda ev,_v=v: _v.get()=="0.000" and _v.set(""))
                e.bind("<FocusOut>", lambda ev,_v=v: _v.set(_v.get() or "0.000"))
                run_entries[ri][s]=v
        sw={"frame":frame,"name":name_var,"runs":run_entries,
            "color":color,"paste_btn":paste_btn,"ref_check":ref_check}
        paste_btn.configure(command=lambda _sw=sw: self._paste_series(_sw))
        self.series_widgets.append(sw)

    def _del_series(self):
        if len(self.series_widgets)<=1: return
        sw=self.series_widgets.pop(); sw["frame"].destroy()

    def _paste_series(self,sw):
        try: text=self.clipboard_get()
        except: messagebox.showinfo("","Pano bos"); return
        rows=parse_paste(text)
        if not rows:
            messagebox.showwarning("","Gecerli veri bulunamadi.\nFormat: 11 sutun"); return
        all_s=["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]
        for ri,row_vals in enumerate(rows[:RUNS_PER_SERIES]):
            for si,val in enumerate(row_vals[:11]):
                s=all_s[si]
                if s in sw["runs"][ri]: sw["runs"][ri][s].set(f"{val:.4f}")

    def _build_right(self,parent):
        self.tabs=ctk.CTkTabview(parent,fg_color="#0e1219",
            segmented_button_fg_color="#1c2336",
            segmented_button_selected_color="#2E75B6",
            segmented_button_unselected_color="#1c2336",
            command=self._on_tab_change)
        self.tabs.pack(fill="both",expand=True,padx=4,pady=4)
        for k in ["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]:
            self.tabs.add(self.T[k])
        self.rf=self.tabs.tab(self.T["tab_results"])
        self.pf=self.tabs.tab(self.T["tab_plot"])
        self.df=self.tabs.tab(self.T["tab_dist"])
        self.sf=self.tabs.tab(self.T["tab_summary"])
        self.cf=self.tabs.tab(self.T["tab_compare"])
        self._last_rendered_tab = None

    def _calculate(self):
        try: self.lo=float(self.e_lo.get()); self.hi=float(self.e_hi.get())
        except: self.lo=15; self.hi=85
        flow=int(self.var_flow.get())
        self.all_series=[]
        for sw in self.series_widgets:
            runs=[]
            for ri in range(RUNS_PER_SERIES):
                m={"Device":0.0}
                for s in DISP_STAGES:
                    try: m[s]=float(sw["runs"][ri][s].get().replace(",","."))
                    except: m[s]=0.0
                r=calc_run(m,flow,self.lo,self.hi,
                    delivered_tp=self.var_delivered_tp.get()); r["run_no"]=ri+1
                runs.append(r)
            avg=calc_series_avg(runs)
            self.all_series.append({
                "name":sw["name"].get(),"color":sw["color"],
                "runs":runs,"avg":avg,
                "is_ref":self.ref_var.get() and (sw==self.series_widgets[0])
            })
        self.lbl_status.configure(text="Hesaplaniyor...")
        self.update_idletasks()
        self._show_results()           # Sonuclar sekmesi hemen render
        self._last_rendered_tab = self.T["tab_results"]
        # Diger sekmeler sekme degisince lazy render edilecek
        # Aktif sekme hangisiyse onu da hemen render et
        cur = self.tabs.get()
        if cur == self.T["tab_plot"]:    self._plot_lp()
        elif cur == self.T["tab_dist"]:  self._plot_dist()
        elif cur == self.T["tab_summary"]: self._show_summary()
        elif cur == self.T["tab_compare"]: self._show_compare()
        self.lbl_status.configure(text=self.T["status_done"])

    def _show_results(self):
        for w in self.rf.winfo_children(): w.destroy()
        scroll=ctk.CTkScrollableFrame(self.rf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        # 5+ seri varsa kumulatif tablo gizli baslar (performans)
        compact_mode = len(self.all_series) >= 5
        for sd in self.all_series:
            rt="  ["+self.T["ref_label"]+"]" if sd["is_ref"] else ""
            ctk.CTkLabel(scroll,text=f"  {sd['name']}{rt}",font=HF,
                text_color=sd["color"],anchor="w").pack(fill="x",padx=8,pady=(10,2))
            for run in sd["runs"]:
                ctk.CTkLabel(scroll,text=f"    Run {run['run_no']}",font=BF,
                    text_color="#aac8e8",anchor="w").pack(fill="x",padx=12,pady=(4,0))
                if "error" in run:
                    ctk.CTkLabel(scroll,text=f"      {self.T['insufficient']} (n={run.get('n',0)})",
                        font=NF,text_color="#ff6060").pack(anchor="w",padx=20); continue
                # Compact modda kumulatif tablo gizli baslar
                tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6)
                if not compact_mode:
                    tf.pack(fill="x",padx=16,pady=2)
                T=self.T; ds=T["dec_sep"]
                hdrs=[T["stage"],"D50",T["mass_mg"],T["cum_mass"],T["cum_pct"],T["valid_pt"],T["probit_z"]]
                ws=[58,66,76,80,68,30,80]
                hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
                for h,w in zip(hdrs,ws):
                    ctk.CTkLabel(hrow,text=h,width=w,font=BF,
                        text_color="white",anchor="center").pack(side="left",padx=1,pady=2)
                vst={v["stage"] for v in run["valid"]}; cum_m=0.0
                for i,row in enumerate(run["cum_data"]):
                    cum_m+=row["mass"]; iv=row["stage"] in vst
                    pz=""
                    if 0<row["u_pct"]<100:
                        try: pz=f"{norm.ppf(row['u_pct']/100):.4f}"
                        except: pass
                    bg="#1a3a1a" if iv else ("#111827" if i%2==0 else "#0e1219")
                    dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
                    _ds=self.T["dec_sep"]
                    vals=[row["stage"],
                          fmt_num(row["d50"],2,_ds) if row["d50"]<900 else "-",
                          fmt_num(row["mass"],4,_ds),
                          fmt_num(cum_m,4,_ds),
                          fmt_num(row["u_pct"],2,_ds),
                          "v" if iv else "",pz]
                    for val,w in zip(vals,ws):
                        ctk.CTkLabel(dr,text=val,width=w,
                            font=BF if iv else NF,
                            text_color="#90ee90" if iv else "#c0d0e0",
                            anchor="center").pack(side="left",padx=1,pady=1)
                pf=ctk.CTkFrame(scroll,fg_color="#1a2540",corner_radius=6)
                pf.pack(fill="x",padx=16,pady=(0,6))
                ds=self.T["dec_sep"]
                params=[
                    (self.T["metered"],    fmt_num(run["metered"],4,ds)),
                    (self.T["delivered"],  fmt_num(run["delivered"],4,ds)),
                    (self.T["fp_dose"],    fmt_num(run["fpd"],4,ds)),
                    (self.T["fp_frac"],    fmt_num(run["fpf"],3,ds)),
                    ("MMAD",               fmt_num(run["mmad"],4,ds)),
                    ("GSD",                fmt_num(run["gsd"],4,ds)),
                    (self.T["slope_lbl"],  fmt_num(run["slope"],4,ds)),
                    (self.T["int_lbl"],    fmt_num(run["intercept"],4,ds)),
                    (self.T["r2_lbl"],     fmt_num(run["r2"],4,ds)),
                    (self.T["n_lbl"],      str(run["n"]))]
                for lbl,val in params:
                    ik=lbl in("FPD","FPF%","MMAD","GSD")
                    ctk.CTkLabel(pf,text=lbl,font=ctk.CTkFont(size=11),
                        text_color="#7090b0").pack(side="left",padx=(8,0))
                    ctk.CTkLabel(pf,text=val,
                        font=ctk.CTkFont(size=12,weight="bold") if ik else ctk.CTkFont(size=12),
                        text_color="#FFC600" if ik else "#e0f0ff").pack(side="left",padx=(2,8))

    def _plot_lp(self):
        for w in self.pf.winfo_children(): w.destroy()
        # Kontrol paneli
        lp_ctrl=ctk.CTkFrame(self.pf,fg_color="#1c2336",corner_radius=6,height=34)
        lp_ctrl.pack(fill="x",padx=4,pady=(4,0)); lp_ctrl.pack_propagate(False)
        ctk.CTkCheckBox(lp_ctrl,text="Sadece Seri Ortalamalari / Series Averages Only",
            variable=self.var_lp_avg_only,font=ctk.CTkFont(size=10),
            text_color="#aac8e8",command=self._plot_lp).pack(side="left",padx=12,pady=6)
        import matplotlib.pyplot as _plt_lp
        _plt_lp.close("all")  # Eski figure'lari temizle
        fig=Figure(figsize=(9,5.2),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        flow=int(self.var_flow.get())
        # 4+ seri varsa otomatik ortalama moda gec (performans)
        avg_only = self.var_lp_avg_only.get() or (len(self.all_series) >= 4)
        for sd in self.all_series:
            if avg_only:
                # Sadece seri ortalaması
                if not sd["avg"]: continue
                valid_runs=[r for r in sd["runs"] if "error" not in r]
                if not valid_runs: continue
                # Ortalama x,y hesapla
                avg_x=np.mean([r["x_reg"] for r in valid_runs],axis=0)
                avg_y=np.mean([r["y_reg"] for r in valid_runs],axis=0)
                lw=3 if sd["is_ref"] else 2
                ax.plot(avg_x,avg_y,"o-",color=sd["color"],alpha=0.9,lw=lw,ms=7,
                    label=sd["name"])
                b=np.sum((avg_x-avg_x.mean())*(avg_y-avg_y.mean()))/np.sum((avg_x-avg_x.mean())**2)
                a=avg_y.mean()-b*avg_x.mean()
                xr=np.linspace(min(avg_x)-0.1,max(avg_x)+0.1,50)
                ax.plot(xr,a+b*xr,"--",color=sd["color"],alpha=0.5,lw=1.5)
            else:
                for run in sd["runs"]:
                    if "error" in run: continue
                    lw=2.5 if sd["is_ref"] else 1.5
                    ax.plot(run["x_reg"],run["y_reg"],"o-",color=sd["color"],
                        alpha=0.85,lw=lw,ms=5,label=f"{sd['name']} R{run['run_no']}")
                    xr=np.linspace(min(run["x_reg"])-0.1,max(run["x_reg"])+0.1,50)
                    ax.plot(xr,run["a"]+run["b"]*xr,"--",color=sd["color"],alpha=0.4,lw=1)
        notes=[]
        for sd in self.all_series:
            if sd["avg"] and "slope" in sd["avg"]["params"]:
                sl=sd["avg"]["params"]["slope"][0]; ic=sd["avg"]["params"]["intercept"][0]
                notes.append(f"{sd['name']}: slope={sl:.3f}  int={ic:.3f}")
        # MMAD dikey çizgisi ekle
        for sd in self.all_series:
            if not sd["avg"] or "mmad" not in sd["avg"]["params"]: continue
            mmad_val = sd["avg"]["params"]["mmad"][0]
            if mmad_val > 0:
                x_mmad = math.log10(mmad_val)
                ax.axvline(x=x_mmad, color=sd["color"], lw=1.2,
                    ls=":", alpha=0.7)
                ax.text(x_mmad, ax.get_ylim()[1]*0.95 if ax.get_ylim()[1]>0 else 0.5,
                    f"MMAD\n{mmad_val:.2f}",
                    color=sd["color"], fontsize=7, ha="center", va="top",
                    bbox=dict(facecolor="#0e1525", alpha=0.5, edgecolor="none", pad=1))
        if notes:
            ax.text(0.02,0.98,"\n".join(notes),transform=ax.transAxes,
                fontsize=9,color="#d0e0f0",va="top",ha="left",
                bbox=dict(facecolor="#0e1525",alpha=0.7,edgecolor="#2a4060",pad=4))
        ax.set_xlabel("log10(D50, um)",color="#7090b0",fontsize=11)
        ax.set_ylabel("Probit z",color="#7090b0",fontsize=11)
        ax.set_title(f"Log-Probit  [{flow} L/min]",color="#FFC600",fontsize=12,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        hdls, lbls = ax.get_legend_handles_labels()
        if hdls:
            n_leg = len(hdls)
            leg_fs = max(6, 9 - max(0, n_leg-6))  # cok seri varsa kucuk font
            ax.legend(fontsize=leg_fs, facecolor="#0e1525", labelcolor="#d0e0f0",
                      loc="best", ncol=2 if n_leg>6 else 1,
                      framealpha=0.85)
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.pf); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)

    def _plot_dist(self):
        for w in self.df.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        # S1 D50=999 olsa bile göster (kütle grafikte görünmeli)
        vis=["S1"]+[s for s in GRAPH_STAGES if s!="S1" and co.get(s,999)<900]
        x=np.arange(len(vis))
        # Limit paneli
        lp=ctk.CTkFrame(self.df,fg_color="#1c2336",corner_radius=6,height=38)
        lp.pack(fill="x",padx=4,pady=(4,0)); lp.pack_propagate(False)
        ctk.CTkLabel(lp,text=self.T["limit_label"],
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#aac8e8").pack(side="left",padx=(8,4),pady=6)
        for val,txt in [("ema",self.T["lim_ema"]),("fda",self.T["lim_fda"]),
                         ("usp",self.T["lim_usp"]),("custom",self.T["lim_custom"])]:
            ctk.CTkRadioButton(lp,text=txt,variable=self.limit_var,value=val,
                font=ctk.CTkFont(size=10),command=self._plot_dist).pack(side="left",padx=5)
        ctk.CTkLabel(lp,text=self.T["lim_pct"],font=ctk.CTkFont(size=10),
            text_color="#aac8e8").pack(side="left",padx=(8,2))
        ctk.CTkEntry(lp,textvariable=self.custom_pct_var,width=46,height=26,
            justify="center",font=ctk.CTkFont(size=11)).pack(side="left",padx=2)
        lim_map={"ema":20,"fda":15,"usp":25}
        try: pct=lim_map.get(self.limit_var.get()) or float(self.custom_pct_var.get())
        except: pct=20
        import matplotlib.pyplot as _plt_dist
        _plt_dist.close("all")
        fig=Figure(figsize=(9,5.0),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        # Throat dahil stage listesi
        vis_all = ["Throat"] + vis
        x_all = np.arange(len(vis_all))
        ref_masses=None; warnings=[]
        for sd in self.all_series:
            if not sd["avg"]: continue
            ms_all=[sd["avg"]["avg_masses"].get(s,0) for s in vis_all]
            valid_runs=[r for r in sd["runs"] if "error" not in r]
            sds_all=[]
            for s in vis_all:
                vals=[r["masses"].get(s,0) for r in valid_runs]
                sds_all.append(float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0)
            lw=3 if sd["is_ref"] else 1.8
            ms_sz=10 if sd["is_ref"] else 5
            ax.plot(x_all,ms_all,color=sd["color"],lw=lw,marker="o",markersize=ms_sz,
                label=sd["name"],zorder=4+(1 if sd["is_ref"] else 0))
            ax.fill_between(x_all,[m-s for m,s in zip(ms_all,sds_all)],
                [m+s for m,s in zip(ms_all,sds_all)],
                color=sd["color"],alpha=0.12,zorder=2)
            ax.errorbar(x_all,ms_all,yerr=sds_all,fmt="none",color=sd["color"],
                capsize=4,lw=1.5,alpha=0.5,zorder=3)
            if sd["is_ref"]: ref_masses=sd["avg"]["avg_masses"]
        if ref_masses:
            rv=[ref_masses.get(s,0) for s in vis_all]
            upper=[v*(1+pct/100) for v in rv]; lower=[v*(1-pct/100) for v in rv]
            ax.plot(x_all,upper,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"+{pct:.0f}%")
            ax.plot(x_all,lower,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"-{pct:.0f}%")
            ax.fill_between(x_all,lower,upper,color="#FF6060",alpha=0.05,zorder=1)
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                mt=[sd["avg"]["avg_masses"].get(s,0) for s in vis_all]
                for s,tv,lo2,hi2 in zip(vis_all,mt,lower,upper):
                    if tv<lo2 or tv>hi2:
                        ds_w=self.T["dec_sep"]
                        yön="yuksek" if tv>hi2 else "dusuk"
                        warnings.append((f"{sd['name']} - {s}: {fmt_num(tv,4,ds_w)} ({yön}) limit", False, False))
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                f2=calc_f2(ref_masses,sd["avg"]["avg_masses"],co)
                if f2 is not None:
                    pf2=self.T["f2_pass"] if f2>=50 else self.T["f2_fail"]
                    ds_f=self.T["dec_sep"]
                    f2_str=fmt_num(f2,1,ds_f)
                    # f2 pass/fail bilgisini tuple olarak sakla: (metin, is_pass)
                    warnings.insert(0,(f"{self.T['f2_label']} {sd['name']}: f2={f2_str} ({pf2})", f2>=50, True))
        ax.set_xticks(list(x_all))
        ax.set_xticklabels(vis_all, rotation=20, ha="right",
            fontsize=10, color="#c0d8f0")
        ax.set_xlabel(self.T["stage"],color="#7090b0",fontsize=11)
        _ylbl = "Ort. Kutle (mg/atis)" if self.lang=="TR" else "Mean Mass (mg/actuation)"
        ax.set_ylabel(_ylbl,color="#7090b0",fontsize=11)
        # Y ekseni virgüllü format
        import matplotlib.ticker as _ticker
        _ds_ax=self.T["dec_sep"]
        ax.yaxis.set_major_formatter(_ticker.FuncFormatter(
            lambda v,p: f"{v:.3f}".replace(".",_ds_ax)))
        if self.lang=="TR":
            ttl=f"APSD  [{flow} L/min]  Ort+/-SD"
            if ref_masses: ttl+=f"  |  Limit +/-{pct:.0f}%"
        else:
            ttl=f"APSD  [{flow} L/min]  Mean+/-SD"
            if ref_masses: ttl+=f"  |  Limit +/-{pct:.0f}%"
        ax.set_title(ttl,color="#FFC600",fontsize=11,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        hdls2, lbls2 = ax.get_legend_handles_labels()
        if hdls2:
            n_leg2 = len(hdls2)
            leg_fs2 = max(6, 9 - max(0, n_leg2-6))
            ax.legend(fontsize=leg_fs2, facecolor="#0e1525", labelcolor="#d0e0f0",
                framealpha=0.85, loc="upper right",
                ncol=2 if n_leg2>6 else 1)
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.df); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)
        if warnings:
            for item in warnings:
                # item = (metin, is_pass, is_f2) veya eski str
                if isinstance(item, tuple):
                    wt, is_pass, is_f2 = item
                else:
                    wt, is_pass, is_f2 = item, False, False
                if is_f2:
                    bg = "#0a2a0a" if is_pass else "#2a0a0a"
                    tc = "#90ee90" if is_pass else "#FF6060"
                else:
                    bg = "#2a0a0a"; tc = "#FFB0B0"
                wf=ctk.CTkFrame(self.df,fg_color=bg,corner_radius=4)
                wf.pack(fill="x",padx=4,pady=1)
                ctk.CTkLabel(wf,text=f"  {wt}",
                    font=ctk.CTkFont(size=11,weight="bold"),
                    text_color=tc,anchor="w").pack(anchor="w",padx=10,pady=3)

    def _show_summary(self):
        for w in self.sf.winfo_children(): w.destroy()
        scroll=ctk.CTkScrollableFrame(self.sf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        try: rsd_lim=float(self.rsd_limit_var.get())
        except: rsd_lim=5.0
        T2=self.T
        params_list=[
            ("metered",T2["metered"]),("delivered",T2["delivered"]),
            ("fpd",T2["fp_dose"]),("fpf",T2["fp_frac"]),
            ("mmad","MMAD (um)"),("gsd","GSD"),
            ("slope",T2["slope_lbl"]),("intercept",T2["int_lbl"]),("r2",T2["r2_lbl"])]
        for sd in self.all_series:
            rt=" ["+self.T["ref_label"]+"]" if sd["is_ref"] else ""
            ctk.CTkLabel(scroll,text=f"  {sd['name']}{rt}",font=HF,
                text_color=sd["color"],anchor="w").pack(fill="x",padx=8,pady=(10,2))
            vr=[r for r in sd["runs"] if "error" not in r]; n=len(vr)
            if n==0:
                ctk.CTkLabel(scroll,text="  Veri yok",font=NF,text_color="#ff6060").pack(anchor="w",padx=20); continue
            tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6); tf.pack(fill="x",padx=12,pady=4)
            cw=[148]+[90]*n+[90,90,82,62]
            T2=self.T; ds=T2["dec_sep"]
            hdrs=[T2["param"]]+[f"{T2['run']} {r['run_no']}" for r in vr]+[T2["mean"],T2["sd"],T2["rsd"],T2["accept"]]
            hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
            for h,w in zip(hdrs,cw):
                ctk.CTkLabel(hrow,text=h,width=w,font=BF,text_color="white",
                    anchor="center").pack(side="left",padx=1,pady=2)
            for i,(key,lbl) in enumerate(params_list):
                vals=[r.get(key) for r in vr if key in r]
                if not vals: continue
                mv=float(np.mean(vals)); sv=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
                rv=sv/mv*100 if mv else 0.0; pf=rv<=rsd_lim
                ik=key in("fpd","fpf","mmad","gsd")
                bg="#1a1a2a" if i%2==0 else "transparent"
                dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
                ctk.CTkLabel(dr,text=lbl,width=cw[0],font=BF if ik else NF,
                    text_color="#FFC600" if ik else "#c0d0e0",anchor="w").pack(side="left",padx=(6,1),pady=2)
                ds=self.T["dec_sep"]
                for r in vr:
                    v=r.get(key,0)
                    ctk.CTkLabel(dr,text=fmt_num(v,4,ds),width=90,font=NF,
                        text_color="#d0e8ff",anchor="center").pack(side="left",padx=1,pady=2)
                for val,w in [(fmt_num(mv,4,ds),90),(fmt_num(sv,4,ds),90),(fmt_num(rv,2,ds),82)]:
                    ctk.CTkLabel(dr,text=val,width=w,font=NF,text_color="#d0e8ff",
                        anchor="center").pack(side="left",padx=1,pady=2)
                ctk.CTkLabel(dr,text="OK" if pf else "FAIL",width=62,font=BF,
                    text_color="#90ee90" if pf else "#ff6060",
                    anchor="center").pack(side="left",padx=1,pady=2)
            dv=[r.get("delivered",0) for r in vr]
            if dv:
                dm=float(np.mean(dv)); ds2=float(np.std(dv,ddof=1)) if len(dv)>1 else 0.0
                dr2=ds2/dm*100 if dm else 0.0
                df=ctk.CTkFrame(scroll,fg_color="#1a2a1a",corner_radius=6)
                df.pack(fill="x",padx=12,pady=(0,4))
                _ds=self.T["dec_sep"]
                ctk.CTkLabel(df,
                    text=f"  {self.T['ddu_label']}: {self.T['ddu_mean']}={fmt_num(dm,4,_ds)}mg  "
                         f"{self.T['ddu_sd']}={fmt_num(ds2,4,_ds)}  "
                         f"{self.T['ddu_rsd']}={fmt_num(dr2,2,_ds)}%  "
                         f"{self.T['cv_label']}={fmt_num(dr2,2,_ds)}%",
                    font=ctk.CTkFont(size=12),text_color="#90ee90").pack(anchor="w",padx=8,pady=4)

    def _show_compare(self):
        for w in self.cf.winfo_children(): w.destroy()
        if not self.all_series: return
        scroll=ctk.CTkScrollableFrame(self.cf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        if len(self.all_series)>=2:
            import matplotlib.pyplot as _plt_cmp; _plt_cmp.close("all")
            fig2=Figure(figsize=(9,3.2),facecolor="#090c12")
            ax1=fig2.add_subplot(121); ax2=fig2.add_subplot(122)
            ax1.set_facecolor("#0e1525"); ax2.set_facecolor("#0e1525")
            names=[sd["name"] for sd in self.all_series if sd["avg"]]
            mmads=[sd["avg"]["params"].get("mmad",(0,))[0] for sd in self.all_series if sd["avg"]]
            gsds=[sd["avg"]["params"].get("gsd",(0,))[0] for sd in self.all_series if sd["avg"]]
            clrs=[sd["color"] for sd in self.all_series if sd["avg"]]
            xi=range(len(names))
            for i,(m,g,c) in enumerate(zip(mmads,gsds,clrs)):
                ax1.bar(i,m,color=c,alpha=0.85,width=0.6)
                ax2.bar(i,g,color=c,alpha=0.85,width=0.6)
            for ax,ttl in [(ax1,self.T["trend_mmad"]),(ax2,self.T["trend_gsd"])]:
                ax.set_xticks(list(xi)); ax.set_xticklabels(names,rotation=20,ha="right",
                    fontsize=9,color="#c0d8f0")
                ax.set_title(ttl,color="#FFC600",fontsize=11,fontweight="bold")
                ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
                ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
            fig2.tight_layout()
            cv2=FigureCanvasTkAgg(fig2,master=scroll); cv2.draw()
            cv2.get_tk_widget().pack(fill="x",pady=(4,0))
        ref_sd=next((sd for sd in self.all_series if sd["is_ref"]),None)
        T3=self.T; ds3=T3["dec_sep"]
        params_list=[
            ("mmad","MMAD (um)"),("gsd","GSD"),
            ("fpd",T3["fp_dose"]),("fpf",T3["fp_frac"]),
            ("slope",T3["slope_lbl"]),("intercept",T3["int_lbl"]),("r2",T3["r2_lbl"])]
        tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6)
        tf.pack(fill="x",padx=8,pady=8)
        cw=[138]+[104]*len(self.all_series)
        hdrs=[self.T["param"]]+[sd["name"]+(" *" if sd["is_ref"] else "") for sd in self.all_series]
        hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
        for h,w in zip(hdrs,cw):
            ctk.CTkLabel(hrow,text=h,width=w,font=BF,text_color="white",
                anchor="center").pack(side="left",padx=1,pady=2)
        for i,(key,lbl) in enumerate(params_list):
            bg="#1a1a2a" if i%2==0 else "transparent"
            dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
            ik=key in("mmad","gsd","fpd","fpf")
            ctk.CTkLabel(dr,text=lbl,width=cw[0],font=BF if ik else NF,
                text_color="#FFC600" if ik else "#c0d0e0",anchor="w").pack(side="left",padx=(6,1),pady=2)
            rv=None
            if ref_sd and ref_sd["avg"] and key in ref_sd["avg"]["params"]:
                rv=ref_sd["avg"]["params"][key][0]
            for sd in self.all_series:
                if not sd["avg"] or key not in sd["avg"]["params"]:
                    ctk.CTkLabel(dr,text="-",width=104,font=NF,text_color="#888",
                        anchor="center").pack(side="left",padx=1,pady=2); continue
                ds3=self.T["dec_sep"]
                val=sd["avg"]["params"][key][0]
                txt=fmt_num(val,4,ds3); clr="#d0e8ff"
                if rv and not sd["is_ref"] and rv>0:
                    diff=(val-rv)/rv*100
                    diff_str=fmt_num(abs(diff),1,ds3)
                    sign="+" if diff>=0 else "-"
                    txt+=f"\n({sign}{diff_str}%)"
                    clr="#90ee90" if abs(diff)<10 else "#FFB060" if abs(diff)<20 else "#FF6060"
                ctk.CTkLabel(dr,text=txt,width=104,font=ctk.CTkFont(size=12),
                    text_color=clr,anchor="center").pack(side="left",padx=1,pady=2)
        if ref_sd and ref_sd["avg"]:
            ctk.CTkLabel(scroll,text=f"  {self.T['f2_label']}",font=HF,
                text_color="#FFD700",anchor="w").pack(fill="x",padx=8,pady=(12,4))
            rm=ref_sd["avg"]["avg_masses"]
            lim_map={"ema":20,"fda":15,"usp":25}
            try: pct=lim_map.get(self.limit_var.get()) or float(self.custom_pct_var.get())
            except: pct=20
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                f2=calc_f2(rm,sd["avg"]["avg_masses"],co)
                if f2 is None: continue
                pf2=self.T["f2_pass"] if f2>=50 else self.T["f2_fail"]
                clr="#90ee90" if f2>=50 else "#FF6060"
                bg_f2="#0a2a0a" if f2>=50 else "#2a0a0a"
                ff=ctk.CTkFrame(scroll,fg_color=bg_f2,corner_radius=6)
                ff.pack(fill="x",padx=12,pady=2)
                ctk.CTkLabel(ff,text=f"  {sd['name']}  f2 = {f2:.1f}   {pf2}",
                    font=ctk.CTkFont(size=13,weight="bold"),text_color=clr).pack(anchor="w",padx=8,pady=6)

    def _toggle_lang(self):
        old_T=self.T; self.lang="EN" if self.lang=="TR" else "TR"; self.T=L[self.lang]
        self.lbl_title.configure(text=self.T["title"])
        self.lbl_sub.configure(text=self.T["subtitle"])
        self.btn_lang.configure(text=self.T["lang_btn"])
        self.btn_add_s.configure(text=self.T["add_series"])
        self.btn_del_s.configure(text=self.T["del_series"])
        self.btn_calc.configure(text=self.T["calculate"])
        self.btn_clr.configure(text=self.T["clear"])
        self.btn_pdf.configure(text=self.T["export_pdf"])
        self.btn_csv.configure(text=self.T["load_csv"])
        self._refresh_cutoffs()
        for sw in self.series_widgets:
            sw["paste_btn"].configure(text=self.T["paste_btn"])
            if sw.get("ref_check"): sw["ref_check"].configure(text=self.T["ref_check"])
        for k in ["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]:
            try:
                self.tabs.rename(old_T[k], self.T[k])
            except Exception:
                pass
        # Aktif sekmeyi yeniden ayarla
        try: self.tabs.set(self.T["tab_results"])
        except: pass

    def _load_csv(self):
        from tkinter import filedialog, messagebox
        path=filedialog.askopenfilename(
            title="CSV Dosyasi Sec / Select CSV File",
            filetypes=[("CSV","*.csv"),("Tum dosyalar / All files","*.*")])
        if not path: return
        try:
            series_dict,flow,warnings=parse_csv(path)
        except Exception as ex:
            messagebox.showerror("CSV Hatasi / Error",str(ex)); return
        if series_dict is None:
            messagebox.showerror("CSV Hatasi",warnings[0] if warnings else "Okunamadi"); return
        # Uyarilar
        if "csv_err_flow" in warnings:
            messagebox.showwarning("",self.T["csv_err_flow"])
        for w in warnings:
            if w.startswith("csv_4runs__"):
                messagebox.showwarning("",self.T["csv_4runs"].format(s=w.replace("csv_4runs__","")))
        # Referans kolonu yoksa popup sor
        has_ref=any(v["ref"] for v in series_dict.values())
        if not has_ref:
            import tkinter as _tk
            names=list(series_dict.keys())
            popup=_tk.Toplevel(self); popup.title(self.T["csv_ref_ask"])
            popup.geometry("440x340"); popup.grab_set()
            popup.configure(bg="#141824")
            ctk.CTkLabel(popup,text=self.T["csv_ref_ask"],
                font=ctk.CTkFont(size=12,weight="bold"),
                text_color="#FFC600").pack(pady=(16,6))
            ref_var=_tk.StringVar(value=self.T["csv_ref_none"])
            sc=ctk.CTkScrollableFrame(popup,fg_color="#1c2336",height=200)
            sc.pack(fill="x",padx=16,pady=4)
            for ch in names+[self.T["csv_ref_none"]]:
                ctk.CTkRadioButton(sc,text=ch,variable=ref_var,value=ch,
                    font=ctk.CTkFont(size=11)).pack(anchor="w",pady=3,padx=8)
            def _ok():
                sel=ref_var.get()
                if sel!=self.T["csv_ref_none"]:
                    for k in series_dict: series_dict[k]["ref"]=(k==sel)
                popup.destroy()
            ctk.CTkButton(popup,text="Tamam / OK",command=_ok,
                fg_color="#1a5a1a",hover_color="#2a8a2a",
                font=ctk.CTkFont(size=11,weight="bold")).pack(pady=10)
            self.wait_window(popup)
        # Flow ayarla
        if flow in NGI_CUTOFFS:
            self.var_flow.set(str(flow)); self._on_flow(str(flow))
        # Mevcut serileri temizle
        for sw in self.series_widgets: sw["frame"].destroy()
        self.series_widgets.clear(); self.all_series=[]
        for tf in [self.rf,self.pf,self.df,self.sf,self.cf]:
            for w in tf.winfo_children(): w.destroy()
        # Serileri yukle
        ref_set=False
        for si,(seri_ad,seri_data) in enumerate(series_dict.items()):
            self._add_series()
            sw=self.series_widgets[-1]
            sw["name"].set(seri_ad)
            # Referans sadece 1. seri icin gecerli
            if si==0:
                self.ref_var.set(seri_data["ref"])
                if seri_data["ref"]: ref_set=True
            # Run degerlerini gir
            for ri,run_data in enumerate(seri_data["runs"][:RUNS_PER_SERIES]):
                for s in DISP_STAGES:
                    val=run_data["masses"].get(s,0.0)
                    sw["runs"][ri][s].set(fmt_num(val,4,"."))
        # 1. seri disinda referans varsa uyar
        ref_list=[k for k,v in series_dict.items() if v["ref"]]
        if ref_list and ref_list[0]!=list(series_dict.keys())[0]:
            messagebox.showinfo("",
                "Referans seri: "+ref_list[0]+
                "\nNot: Referans secenegi yalnizca 1. seride etkin.\n"
                "Referans seriyi 1. siraya tasiyin.")
        n_s=len(series_dict)
        n_r=sum(len(v["runs"]) for v in series_dict.values())
        self.lbl_status.configure(text=self.T["csv_loaded"].format(n=n_s,r=n_r))

    def _clear(self):
        for sw in self.series_widgets:
            for ri in range(RUNS_PER_SERIES):
                for s in sw["runs"][ri]: sw["runs"][ri][s].set("0.000")
        self.all_series=[]
        for tf in [self.rf,self.pf,self.df,self.sf,self.cf]:
            for w in tf.winfo_children(): w.destroy()

    def _export_pdf(self):
        if not self.all_series:
            messagebox.showwarning("","Oncelikle hesaplama yapiniz."); return
        path=filedialog.asksaveasfilename(defaultextension=".pdf",
            filetypes=[("PDF","*.pdf")],
            initialfile=f"NGI_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf")
        if not path: return
        meta={"product":self.e_product.get(),"batch":self.e_batch.get(),
              "operator":self.e_operator.get(),"date":self.e_date.get()}
        lm={"ema":20,"fda":15,"usp":25}
        try: pct=lm.get(self.limit_var.get()) or float(self.custom_pct_var.get())
        except: pct=20
        try:
            make_pdf_multi(path,self.all_series,meta,int(self.var_flow.get()),self.T,pct,
                           rsd_lim=float(self.rsd_limit_var.get()) if self.rsd_limit_var.get() else 5.0,
                           lang=self.lang)
            # Otomatik aç
            import subprocess, platform
            try:
                if platform.system()=="Windows": os.startfile(path)
                elif platform.system()=="Darwin": subprocess.Popen(["open",path])
                else: subprocess.Popen(["xdg-open",path])
            except: pass
            messagebox.showinfo("",f"PDF kaydedildi:\n{path}")
        except Exception as ex:
            import traceback as _tb
            full = _tb.format_exc()
            messagebox.showerror("PDF Hatasi",
                str(ex) + "\n\n--- Detay ---\n" + full[-500:])

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
                 ("fpd","FPD [mg]"),("fpf","FPF [%]"),
                 ("mmad","MMAD [um]"),("gsd","GSD"),
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
            is_key = k in ("fpd","fpf","mmad","gsd")
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
    BW_LS     = ["-","--","-.",":",(0,(5,1)),(0,(3,1,1,1))]

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
        ax_ap.plot(x_ap, ms_ap, "o"+ls, color=col, lw=lw, ms=5,
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


if __name__=="__main__":
    app=NGIApp(); app.mainloop()
