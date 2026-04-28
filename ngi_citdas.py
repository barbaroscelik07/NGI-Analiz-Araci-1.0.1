"""NGI Cascade Impactor Analysis Tool v4 - Ph.Eur 2.9.18 / USP <601>
Çoklu seri desteği, her seri 3 run, Excel yapıştırma, karşılaştırma grafiği
"""
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
    base = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(
        sys.executable if getattr(sys,'frozen',False) else __file__)))
    return os.path.join(base, rel)

# ── NGI Cut-off D50 (µm) ──────────────────────────────────────────────────────
NGI_CUTOFFS = {
    15:  {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":8.61,"S3":5.39,"S4":3.30,"S5":2.08,"S6":1.36,"S7":0.98,"MOC":0.54},
    30:  {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":4.46,"S3":2.82,"S4":1.66,"S5":0.94,"S6":0.55,"S7":0.34,"MOC":0.21},
    45:  {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":3.42,"S3":2.09,"S4":1.21,"S5":0.72,"S6":0.40,"S7":0.24,"MOC":0.13},
    60:  {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":8.06,"S3":4.46,"S4":2.82,"S5":1.66,"S6":0.94,"S7":0.55,"MOC":0.34},
    90:  {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":2.08,"S3":1.36,"S4":0.98,"S5":0.55,"S6":0.34,"S7":0.21,"MOC":0.10},
    100: {"Device":999,"Throat":999,"Presep":999,"S1":999,
          "S2":1.78,"S3":1.12,"S4":0.69,"S5":0.43,"S6":0.25,"S7":0.14,"MOC":0.08},
}
STAGE_ORDER = ["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7"]
ALL_KEYS    = STAGE_ORDER + ["MOC"]
ISM_STAGES  = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
RUNS_PER_SERIES = 3
CP = ["#2E75B6","#ED7D31","#70AD47","#E84040","#7030A0",
      "#00B0F0","#D4A000","#C00000","#00B050","#FF69B4"]

L = {
"TR":{"title":"NGI \u0130mpaktor Analiz Arac\u0131","subtitle":"Ph.Eur 2.9.18 / USP <601>",
 "lang_btn":"English","product":"\u00dcr\u00fcn Ad\u0131","batch":"Lot No.",
 "operator":"Analist","date":"Tarih","flow_rate":"Ak\u0131\u015f H\u0131z\u0131",
 "add_series":"+ Seri Ekle","del_series":"Seri Sil",
 "calculate":"Hesapla","clear":"Temizle","export_pdf":"PDF Rapor",
 "tab_results":"Sonu\u00e7lar","tab_plot":"Log-Probit",
 "tab_dist":"Da\u011f\u0131l\u0131m","tab_summary":"\u00d6zet","tab_compare":"Kar\u015f\u0131la\u015ft\u0131rma",
 "series":"Seri","run":"Run","series_name":"Seri Ad\u0131",
 "paste_hint":"Excel'den kopyalay\u0131p yap\u0131\u015ft\u0131r (Ctrl+V)",
 "paste_btn":"Yapıştır","mmad":"MMAD (\u00b5m)","gsd":"GSD",
 "fpd":"FPD (mg)","fpf":"FPF (%)","r2":"R\u00b2","n":"n",
 "slope":"Slope","intercept":"Intercept",
 "mean":"Ort.","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Yetersiz nokta","status_ready":"Haz\u0131r.",
 "status_done":"Hesaplama tamamland\u0131.",
 "err_val":"Ge\u00e7ersiz de\u011fer","err_nodata":"Veri yok",
 "valid_range":"Ge\u00e7erlilik (%)",
 "cutoff_title":"Cut-off D50 (\u00b5m)"},
"EN":{"title":"NGI Cascade Impactor Analysis","subtitle":"Ph.Eur 2.9.18 / USP <601>",
 "lang_btn":"T\u00fcrk\u00e7e","product":"Product","batch":"Batch No.",
 "operator":"Analyst","date":"Date","flow_rate":"Flow Rate",
 "add_series":"+ Add Series","del_series":"Del Series",
 "calculate":"Calculate","clear":"Clear","export_pdf":"PDF Report",
 "tab_results":"Results","tab_plot":"Log-Probit",
 "tab_dist":"Distribution","tab_summary":"Summary","tab_compare":"Compare",
 "series":"Series","run":"Run","series_name":"Series Name",
 "paste_hint":"Copy from Excel and paste (Ctrl+V)",
 "paste_btn":"Paste","mmad":"MMAD (\u00b5m)","gsd":"GSD",
 "fpd":"FPD (mg)","fpf":"FPF (%)","r2":"R\u00b2","n":"n",
 "slope":"Slope","intercept":"Intercept",
 "mean":"Mean","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Insufficient pts","status_ready":"Ready.",
 "status_done":"Calculation complete.",
 "err_val":"Invalid value","err_nodata":"No data",
 "valid_range":"Valid Range (%)",
 "cutoff_title":"Cut-off D50 (\u00b5m)"},
}

# ── Hesaplama motoru ───────────────────────────────────────────────────────────
def calc_run(masses, flow, lo=15, hi=85):
    co  = NGI_CUTOFFS[flow]
    ism = sum(masses.get(s,0) for s in ISM_STAGES)
    metered = sum(masses.get(s,0) for s in ALL_KEYS)
    if ism <= 0:
        return {"error":"no_data","metered":metered,"delivered":ism,"masses":masses}
    cum=[]
    for s in ALL_KEYS:
        u=0.0
        if s in ISM_STAGES:
            u = sum(masses.get(x,0) for x in ISM_STAGES
                    if co.get(x,999) <= co.get(s,999)) / ism * 100
        cum.append({"stage":s,"d50":co.get(s,999),"mass":masses.get(s,0),"u_pct":u})
    valid=[r for r in cum if r["stage"] in ISM_STAGES
           and co.get(r["stage"],999)<900 and lo<r["u_pct"]<hi]
    res={"metered":metered,"delivered":ism,"cum_data":cum,
         "valid":valid,"masses":masses,"flow":flow}
    if len(valid)<2:
        res["error"]="insufficient"; res["n"]=len(valid); return res
    x=np.array([math.log10(v["d50"]) for v in valid])
    y=np.array([norm.ppf(v["u_pct"]/100) for v in valid])
    b=np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
    a=y.mean()-b*x.mean()
    yp=a+b*x; ss_r=np.sum((y-yp)**2); ss_t=np.sum((y-y.mean())**2)
    r2=1-ss_r/ss_t if ss_t>0 else 1.0
    z5=a+b*math.log10(5)
    res.update({"n":len(valid),"a":a,"b":b,"slope":b,"intercept":a+5,"r2":r2,
                "mmad":10**(-a/b),"gsd":10**(1/b),
                "fpd":norm.cdf(z5)*ism,"fpf":norm.cdf(z5)*ism/metered*100,
                "x_reg":x,"y_reg":y})
    return res

def calc_series_avg(runs):
    """3 run'ın ortalama kütle ve parametrelerini hesapla"""
    valid=[r for r in runs if "error" not in r]
    if not valid: return None
    # Ortalama kütleler (stage bazında)
    avg_masses={}
    for s in ALL_KEYS:
        vals=[r["masses"].get(s,0) for r in valid]
        avg_masses[s]=float(np.mean(vals))
    # Ortalama parametreler
    params={}
    for p in ["mmad","gsd","fpd","fpf","metered","delivered","slope","intercept","r2"]:
        vals=[r[p] for r in valid if p in r]
        if vals:
            m=float(np.mean(vals)); s=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            params[p]=(m,s,s/m*100 if m else 0.0)
    return {"avg_masses":avg_masses,"params":params,"n_valid":len(valid)}

def parse_paste(text):
    """Excel'den yapıştırılan metni parse et — başlıklı veya başlıksız"""
    lines=[l.strip() for l in text.strip().splitlines() if l.strip()]
    if not lines: return None
    # Başlık satırı var mı? (sayısal olmayan ilk token)
    def is_header(line):
        first=line.split('\t')[0].strip()
        try: float(first.replace(',','.')); return False
        except: return True
    if is_header(lines[0]): lines=lines[1:]
    if not lines: return None
    result=[]
    for line in lines:
        tokens=[t.strip() for t in line.split('\t')]
        try:
            vals=[float(t.replace(',','.')) for t in tokens if t]
            if len(vals)>=11:
                result.append(vals[:11])  # Device,Throat,Presep,S1-S7,MOC
        except: pass
    return result if result else None

# ── Ana uygulama ──────────────────────────────────────────────────────────────
class NGIApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.lang="TR"; self.T=L["TR"]
        self.all_series=[]   # Her eleman: {"name":str, "runs":[3 run result], "avg":dict}
        self.flow=30; self.lo=15; self.hi=85
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")
        self.title(self.T["title"])
        self.geometry("1400x920"); self.minsize(1100,750)
        ico=resource_path("icon.ico")
        if os.path.exists(ico):
            try: self.iconbitmap(ico)
            except: pass
        self._build_ui()
        self.lbl_status.configure(text=self.T["status_ready"])

    # ── UI inşa ────────────────────────────────────────────────────────────────
    def _build_ui(self):
        # Header
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
        # Body
        body=ctk.CTkFrame(self,fg_color="transparent"); body.pack(fill="both",expand=True)
        # Sol panel
        self.left=ctk.CTkScrollableFrame(body,width=440,fg_color="#141824",corner_radius=0)
        self.left.pack(side="left",fill="y")
        self._build_left()
        # Sağ panel
        right=ctk.CTkFrame(body,fg_color="#0e1219"); right.pack(side="left",fill="both",expand=True)
        self._build_right(right)
        # Status bar
        sb=ctk.CTkFrame(self,height=24,fg_color="#090c12",corner_radius=0)
        sb.pack(fill="x",side="bottom"); sb.pack_propagate(False)
        self.lbl_status=ctk.CTkLabel(sb,text="",anchor="w",
            font=ctk.CTkFont(size=11),text_color="#4a7a9a")
        self.lbl_status.pack(side="left",padx=10)

    def _lbl(self,text,color="#7a9abf"):
        ctk.CTkLabel(self.left,text=text,anchor="w",
            font=ctk.CTkFont(size=11,weight="bold"),text_color=color).pack(
            fill="x",padx=8,pady=(10,2))

    def _card(self):
        f=ctk.CTkFrame(self.left,fg_color="#1c2336",corner_radius=8)
        f.pack(fill="x",padx=6,pady=3)
        i=ctk.CTkFrame(f,fg_color="transparent"); i.pack(fill="x",padx=10,pady=8)
        return i

    def _build_left(self):
        # Meta
        self._lbl("● Bilgi","#FFC600")
        mc=self._card(); self.mv={}
        for k,d in [("product",""),("batch",""),("operator",""),
                    ("date",datetime.now().strftime("%d.%m.%Y"))]:
            r=ctk.CTkFrame(mc,fg_color="transparent"); r.pack(fill="x",pady=2)
            ctk.CTkLabel(r,text=self.T.get(k,k),width=80,anchor="w",
                font=ctk.CTkFont(size=11),text_color="#7a9abf").pack(side="left")
            v=ctk.StringVar(value=d); self.mv[k]=v
            ctk.CTkEntry(r,textvariable=v,height=26,font=ctk.CTkFont(size=11)
                ).pack(side="left",fill="x",expand=True,padx=4)
        # Flow + ayarlar
        self._lbl("● "+self.T["flow_rate"],"#FFC600")
        fc=self._card()
        fr=ctk.CTkFrame(fc,fg_color="transparent"); fr.pack(fill="x",pady=2)
        ctk.CTkLabel(fr,text=self.T["flow_rate"],width=100,anchor="w",
            font=ctk.CTkFont(size=11),text_color="#7a9abf").pack(side="left")
        self.var_flow=ctk.StringVar(value="30")
        ctk.CTkOptionMenu(fr,values=["15","30","45","60","90","100"],
            variable=self.var_flow,command=self._on_flow,width=90,height=26,
            font=ctk.CTkFont(size=12,weight="bold")).pack(side="left",padx=6)
        ctk.CTkLabel(fr,text="L/min",font=ctk.CTkFont(size=11),
            text_color="#7a9abf").pack(side="left")
        self.cbox=ctk.CTkFrame(fc,fg_color="#090c14",corner_radius=6)
        self.cbox.pack(fill="x",pady=(6,2)); self._refresh_cutoffs()
        # Geçerlilik aralığı
        r2=ctk.CTkFrame(fc,fg_color="transparent"); r2.pack(fill="x",pady=4)
        self.lbl_vr=ctk.CTkLabel(r2,text=self.T["valid_range"],width=130,
            anchor="w",font=ctk.CTkFont(size=11),text_color="#7a9abf")
        self.lbl_vr.pack(side="left")
        self.var_lo=ctk.StringVar(value="15"); self.var_hi=ctk.StringVar(value="85")
        ctk.CTkEntry(r2,textvariable=self.var_lo,width=48,height=26).pack(side="left",padx=2)
        ctk.CTkLabel(r2,text="–").pack(side="left",padx=2)
        ctk.CTkEntry(r2,textvariable=self.var_hi,width=48,height=26).pack(side="left",padx=2)
        # Seriler
        self._lbl("● Seriler (her seri = 3 Run)","#FFC600")
        self.series_box=ctk.CTkFrame(self.left,fg_color="transparent")
        self.series_box.pack(fill="x",padx=4,pady=2)
        self.series_widgets=[]
        self._add_series()
        # Butonlar
        br=ctk.CTkFrame(self.left,fg_color="transparent"); br.pack(fill="x",padx=6,pady=4)
        self.btn_add_s=ctk.CTkButton(br,text=self.T["add_series"],command=self._add_series,
            height=30,width=115,fg_color="#1a4d30",hover_color="#256640")
        self.btn_add_s.pack(side="left",padx=3)
        self.btn_del_s=ctk.CTkButton(br,text=self.T["del_series"],command=self._del_series,
            height=30,width=100,fg_color="#4d1a1a",hover_color="#662525")
        self.btn_del_s.pack(side="left",padx=3)
        ar=ctk.CTkFrame(self.left,fg_color="transparent"); ar.pack(fill="x",padx=6,pady=8)
        self.btn_calc=ctk.CTkButton(ar,text=self.T["calculate"],command=self._calculate,
            height=38,width=140,font=ctk.CTkFont(size=13,weight="bold"),
            fg_color="#002D62",hover_color="#003a80")
        self.btn_calc.pack(side="left",padx=3)
        self.btn_clr=ctk.CTkButton(ar,text=self.T["clear"],command=self._clear,
            height=38,width=90,fg_color="#333",hover_color="#555")
        self.btn_clr.pack(side="left",padx=3)
        self.btn_pdf=ctk.CTkButton(ar,text=self.T["export_pdf"],command=self._export_pdf,
            height=38,width=120,fg_color="#1a4d30",hover_color="#256640")
        self.btn_pdf.pack(side="left",padx=3)

    def _build_right(self,parent):
        self.tabs=ctk.CTkTabview(parent,fg_color="#0e1219",
            segmented_button_fg_color="#141824",
            segmented_button_selected_color="#002D62",
            segmented_button_unselected_color="#1c2336")
        self.tabs.pack(fill="both",expand=True,padx=4,pady=4)
        for k in ["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]:
            self.tabs.add(self.T[k])
        self.tabs.set(self.T["tab_results"])
        self.rf=ctk.CTkScrollableFrame(self.tabs.tab(self.T["tab_results"]),fg_color="transparent")
        self.rf.pack(fill="both",expand=True)
        self.pf=ctk.CTkFrame(self.tabs.tab(self.T["tab_plot"]),fg_color="#090c12")
        self.pf.pack(fill="both",expand=True)
        self.df=ctk.CTkFrame(self.tabs.tab(self.T["tab_dist"]),fg_color="#090c12")
        self.df.pack(fill="both",expand=True)
        self.sf=ctk.CTkScrollableFrame(self.tabs.tab(self.T["tab_summary"]),fg_color="transparent")
        self.sf.pack(fill="both",expand=True)
        self.cf=ctk.CTkFrame(self.tabs.tab(self.T["tab_compare"]),fg_color="#090c12")
        self.cf.pack(fill="both",expand=True)

    # ── Cut-off paneli ─────────────────────────────────────────────────────────
    def _refresh_cutoffs(self):
        for w in self.cbox.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        ctk.CTkLabel(self.cbox,text=f" {self.T['cutoff_title']}  [{flow} L/min]",
            font=ctk.CTkFont(size=10,weight="bold"),text_color="#FFC600",anchor="w"
            ).pack(fill="x",padx=6,pady=(4,2))
        gf=ctk.CTkFrame(self.cbox,fg_color="transparent"); gf.pack(padx=6,pady=(0,6))
        visible=[s for s in ALL_KEYS if co.get(s,999)<900]
        for s in visible:
            c=ctk.CTkFrame(gf,fg_color="#001a40",corner_radius=4,width=46)
            c.pack(side="left",padx=2); c.pack_propagate(False)
            ctk.CTkLabel(c,text=s,font=ctk.CTkFont(size=8),text_color="#aac8e8").pack(pady=(2,0))
            ctk.CTkLabel(c,text=f"{co[s]:.2f}",font=ctk.CTkFont(size=9,weight="bold"),
                text_color="#FFC600").pack(pady=(0,2))

    def _on_flow(self,v): self.flow=int(v); self._refresh_cutoffs()

    # ── Seri widget'ı ──────────────────────────────────────────────────────────
    def _add_series(self):
        idx=len(self.series_widgets)+1
        color=CP[(idx-1)%len(CP)]
        frame=ctk.CTkFrame(self.series_box,fg_color="#1c2336",corner_radius=8)
        frame.pack(fill="x",pady=4,padx=2)
        # Seri başlık + isim
        hf=ctk.CTkFrame(frame,fg_color="#001a40",corner_radius=6)
        hf.pack(fill="x",padx=6,pady=(6,4))
        ctk.CTkLabel(hf,text=f"  {self.T['series']} {idx}",
            font=ctk.CTkFont(size=11,weight="bold"),text_color=color
            ).pack(side="left",pady=4)
        name_var=ctk.StringVar(value=f"Seri {idx}")
        ctk.CTkEntry(hf,textvariable=name_var,height=24,width=120,
            font=ctk.CTkFont(size=10),placeholder_text="Seri adi"
            ).pack(side="left",padx=8)
        paste_btn=ctk.CTkButton(hf,text=self.T["paste_btn"],width=80,height=24,
            font=ctk.CTkFont(size=10),fg_color="#003580",hover_color="#0055c0")
        paste_btn.pack(side="right",padx=6)

        # Grid: satır=stage, sütun=run
        # Üst satır: boş + Run1, Run2, Run3 başlıkları
        grid=ctk.CTkFrame(frame,fg_color="transparent")
        grid.pack(fill="x",padx=8,pady=(4,6))

        # Başlık satırı
        ctk.CTkLabel(grid,text="Stage",width=58,
            font=ctk.CTkFont(size=9,weight="bold"),
            text_color="#5a8ab0",anchor="w").grid(row=0,column=0,padx=2,pady=1)
        for ri in range(RUNS_PER_SERIES):
            ctk.CTkLabel(grid,text=f"Run {ri+1}",width=72,
                font=ctk.CTkFont(size=9,weight="bold"),
                text_color=color,anchor="center").grid(row=0,column=ri+1,padx=2,pady=1)

        # Stage satırları
        run_entries=[{} for _ in range(RUNS_PER_SERIES)]
        for si,s in enumerate(ALL_KEYS):
            row_i=si+1
            # Stage etiketi
            ctk.CTkLabel(grid,text=s,width=58,
                font=ctk.CTkFont(size=9),text_color="#aac8e8",anchor="w"
                ).grid(row=row_i,column=0,padx=2,pady=1)
            # Her run için giriş kutusu
            for ri in range(RUNS_PER_SERIES):
                v=ctk.StringVar(value="0.000")
                e=ctk.CTkEntry(grid,textvariable=v,height=22,width=72,
                    font=ctk.CTkFont(size=10),justify="center")
                e.grid(row=row_i,column=ri+1,padx=2,pady=1)
                e.bind("<FocusIn>", lambda ev,_v=v: _v.get()=="0.000" and _v.set(""))
                e.bind("<FocusOut>",lambda ev,_v=v: _v.set(_v.get() or "0.000"))
                run_entries[ri][s]=v

        sw={"frame":frame,"name":name_var,"runs":run_entries,
            "color":color,"paste_btn":paste_btn}
        paste_btn.configure(command=lambda _sw=sw: self._paste_series(_sw))
        self.series_widgets.append(sw)

    def _del_series(self):
        if len(self.series_widgets)<=1: return
        self.series_widgets.pop()["frame"].destroy()

    def _paste_series(self,sw):
        """Clipboard'dan veri yapıştır"""
        try:
            text=self.clipboard_get()
        except:
            messagebox.showinfo("","Pano bo\u015f"); return
        rows=parse_paste(text)
        if not rows:
            messagebox.showwarning("",
                "Ge\u00e7erli veri bulunamad\u0131.\n"
                "Format: 11 s\u00fctun (Device,Throat,Presep,S1-S7,MOC)\n"
                "Ba\u015fl\u0131kl\u0131 veya ba\u015fl\u0131ks\u0131z olabilir.")
            return
        for ri,row_vals in enumerate(rows[:RUNS_PER_SERIES]):
            for si,s in enumerate(ALL_KEYS):
                if si<len(row_vals):
                    sw["runs"][ri][s].set(f"{row_vals[si]:.4f}")
        self.lbl_status.configure(
            text=f"{sw['name'].get()}: {min(len(rows),RUNS_PER_SERIES)} run yap\u0131\u015ft\u0131r\u0131ld\u0131.")

    # ── Hesaplama ──────────────────────────────────────────────────────────────
    def _calculate(self):
        try: lo=float(self.var_lo.get()); hi=float(self.var_hi.get())
        except: lo,hi=15,85
        flow=int(self.var_flow.get())
        self.all_series=[]
        for si,sw in enumerate(self.series_widgets):
            runs=[]
            for ri,run_vars in enumerate(sw["runs"]):
                try:
                    m={s:float(run_vars[s].get().replace(",",".")) for s in ALL_KEYS}
                except:
                    messagebox.showerror("",
                        f"{self.T['err_val']} — {sw['name'].get()} Run {ri+1}")
                    return
                r=calc_run(m,flow,lo,hi)
                r["run_no"]=ri+1; r["series_no"]=si+1
                runs.append(r)
            avg=calc_series_avg(runs)
            self.all_series.append({
                "name":sw["name"].get(),"runs":runs,"avg":avg,
                "color":sw["color"],"flow":flow
            })
        self._show_results()
        self.after(100,self._plot_lp)
        self.after(200,self._plot_dist)
        self.after(300,self._show_summary)
        self.after(400,self._plot_compare)
        self.tabs.set(self.T["tab_results"])
        self.lbl_status.configure(text=self.T["status_done"])

    # ── Sonuçlar sekmesi ───────────────────────────────────────────────────────
    def _show_results(self):
        for w in self.rf.winfo_children(): w.destroy()
        for sd in self.all_series:
            c=sd["color"]
            # Seri başlık
            hf=ctk.CTkFrame(self.rf,fg_color=c,corner_radius=6,height=28)
            hf.pack(fill="x",pady=(8,2),padx=6); hf.pack_propagate(False)
            ctk.CTkLabel(hf,text=f"  {sd['name']}  —  {sd['flow']} L/min",
                font=ctk.CTkFont(size=13,weight="bold"),text_color="white",anchor="w"
                ).pack(side="left",padx=10)
            # 3 run yan yana
            runs_frame=ctk.CTkFrame(self.rf,fg_color="#1c2336",corner_radius=8)
            runs_frame.pack(fill="x",padx=6,pady=3)
            for run in sd["runs"]:
                rf2=ctk.CTkFrame(runs_frame,fg_color="#0f1628",corner_radius=6)
                rf2.pack(side="left",fill="both",expand=True,padx=4,pady=6)
                ctk.CTkLabel(rf2,text=f"Run {run['run_no']}",
                    font=ctk.CTkFont(size=11,weight="bold"),text_color=c
                    ).pack(pady=(6,2))
                if "error" in run:
                    ctk.CTkLabel(rf2,text=self.T["insufficient"],
                        font=ctk.CTkFont(size=10),text_color="#ff6644").pack()
                else:
                    params=[
                        ("Metered",f"{run['metered']:.4f}",False),
                        ("Delivered",f"{run['delivered']:.4f}",False),
                        ("FPD",f"{run['fpd']:.4f}",True),
                        ("FPF%",f"{run['fpf']:.3f}",True),
                        ("MMAD",f"{run['mmad']:.4f}",True),
                        ("GSD",f"{run['gsd']:.4f}",True),
                        ("Slope",f"{run['slope']:.4f}",False),
                        ("Intercept",f"{run['intercept']:.4f}",False),
                        ("R²",f"{run['r2']:.4f}",False),
                        ("n",str(run['n']),False),
                    ]
                    for lbl,val,hl in params:
                        pr=ctk.CTkFrame(rf2,fg_color="transparent")
                        pr.pack(fill="x",padx=6,pady=1)
                        ctk.CTkLabel(pr,text=lbl,width=70,anchor="w",
                            font=ctk.CTkFont(size=10),text_color="#7a9abf"
                            ).pack(side="left")
                        ctk.CTkLabel(pr,text=val,anchor="e",
                            font=ctk.CTkFont(size=10,
                                weight="bold" if hl else "normal"),
                            text_color="#FFC600" if hl else "#e0e0e0"
                            ).pack(side="right",padx=4)
                ctk.CTkFrame(rf2,height=4,fg_color="transparent").pack()
            # Seri ortalaması
            if sd["avg"]:
                af=ctk.CTkFrame(self.rf,fg_color="#0a1422",corner_radius=6)
                af.pack(fill="x",padx=6,pady=(0,6))
                ctk.CTkLabel(af,text=f"  {sd['name']} Ortalama",
                    font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFC600"
                    ).pack(anchor="w",padx=10,pady=(6,2))
                pr2=sd["avg"]["params"]
                row=ctk.CTkFrame(af,fg_color="transparent"); row.pack(fill="x",padx=10,pady=(0,6))
                for key,lbl in [("fpd","FPD"),("fpf","FPF%"),("mmad","MMAD"),("gsd","GSD")]:
                    if key in pr2:
                        m,s,rsd=pr2[key]
                        col=ctk.CTkFrame(row,fg_color="#1c2336",corner_radius=6)
                        col.pack(side="left",padx=3,pady=2,fill="x",expand=True)
                        ctk.CTkLabel(col,text=lbl,font=ctk.CTkFont(size=9),
                            text_color="#7a9abf").pack(pady=(4,0))
                        ctk.CTkLabel(col,text=f"{m:.4f}",
                            font=ctk.CTkFont(size=11,weight="bold"),
                            text_color="#FFC600").pack()
                        ctk.CTkLabel(col,text=f"±{s:.4f}",
                            font=ctk.CTkFont(size=9),text_color="#aaa").pack(pady=(0,4))

    # ── Log-Probit grafiği ─────────────────────────────────────────────────────
    def _plot_lp(self):
        for w in self.pf.winfo_children(): w.destroy()
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        markers=['o','s','^','D','v','<']
        for sd in self.all_series:
            c=sd["color"]
            for run in sd["runs"]:
                if "error" in run: continue
                lbl=f"{sd['name']} R{run['run_no']}"
                ax.scatter(run["x_reg"],run["y_reg"],color=c,s=50,zorder=5,
                    marker=markers[run["run_no"]-1],alpha=0.7)
                xl=np.linspace(min(run["x_reg"])-0.1,max(run["x_reg"])+0.1,100)
                ax.plot(xl,run["a"]+run["b"]*xl,color=c,lw=1,ls="--",alpha=0.5)
            # Seri ortalaması kalın çizgi
            if sd["avg"]:
                valid_runs=[r for r in sd["runs"] if "error" not in r]
                if len(valid_runs)>=2:
                    xs=np.concatenate([r["x_reg"] for r in valid_runs])
                    ys=np.concatenate([r["y_reg"] for r in valid_runs])
                    xr=np.linspace(xs.min()-0.1,xs.max()+0.1,100)
                    b=np.sum((xs-xs.mean())*(ys-ys.mean()))/np.sum((xs-xs.mean())**2)
                    a=ys.mean()-b*xs.mean()
                    ax.plot(xr,a+b*xr,color=c,lw=2.5,label=sd["name"])
                    mx=math.log10(sd["avg"]["params"]["mmad"][0])
                    ax.axvline(mx,color=c,ls=":",lw=1.5,alpha=0.8)
                    ax.text(mx+0.01,ax.get_ylim()[0] if ax.get_ylim()[0]>-3.5 else -3.2,
                        f"MMAD={sd['avg']['params']['mmad'][0]:.3f}",
                        color=c,fontsize=7.5,rotation=90,va="bottom")
        ax.set_xlabel("log₁₀(D50, µm)",color="#7090b0",fontsize=10)
        ax.set_ylabel("Probit z",color="#7090b0",fontsize=10)
        ax.set_title("Log-Probit Grafiği",color="#FFC600",fontsize=12,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        # Slope ve intercept annotasyonu
        note_lines=[]
        for sd in self.all_series:
            if sd["avg"] and "slope" in sd["avg"]["params"]:
                sl=sd["avg"]["params"]["slope"][0]
                ic=sd["avg"]["params"]["intercept"][0]
                note_lines.append(f"{sd['name']}: slope={sl:.3f}  intercept={ic:.3f}")
        if note_lines:
            ax.text(0.02,0.98,"\n".join(note_lines),transform=ax.transAxes,
                fontsize=8,color="#d0e0f0",va="top",ha="left",
                bbox=dict(facecolor="#0e1525",alpha=0.7,edgecolor="#2a4060",pad=4))
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.pf); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)

    # ── Dağılım grafiği ────────────────────────────────────────────────────────
    def _plot_dist(self):
        for w in self.df.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        # Sadece gecerli D50 olan stagelar (APSD grafigi gibi)
        vis_stages=[s for s in ALL_KEYS if co.get(s,999)<900]
        x=np.arange(len(vis_stages))
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        for sd in self.all_series:
            if not sd["avg"]: continue
            valid_runs=[r for r in sd["runs"] if "error" not in r]
            ms=[sd["avg"]["avg_masses"].get(s,0) for s in vis_stages]
            sds=[]
            for s in vis_stages:
                vals=[r["masses"].get(s,0) for r in valid_runs]
                sds.append(float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0)
            # APSD tarzı: çizgi + işaretçi + hata çubuğu
            ax.plot(x, ms, color=sd["color"], lw=2, marker="o",
                markersize=6, label=sd["name"], zorder=4)
            ax.fill_between(x,
                [m-s for m,s in zip(ms,sds)],
                [m+s for m,s in zip(ms,sds)],
                color=sd["color"], alpha=0.15, zorder=2)
            ax.errorbar(x, ms, yerr=sds, fmt="none",
                color=sd["color"], capsize=4, lw=1.5, alpha=0.6, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(vis_stages, rotation=0, fontsize=10, color="#c0d8f0")
        ax.set_xlabel("Stage", color="#7090b0", fontsize=11)
        ax.set_ylabel("Ortalama Kütle (mg/atış)",
            color="#7090b0", fontsize=11)
        ax.set_title(
            f"APSD — Stage Kütle Dağılımı  [{flow} L/min]"
            f"  |  Ortalama ± SD",
            color="#FFC600", fontsize=12, fontweight="bold")
        ax.tick_params(colors="#7090b0")
        ax.spines[:].set_color("#2a4060")
        ax.grid(True, color="#1a3050", ls="--", alpha=0.5)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=10, facecolor="#0e1525", labelcolor="#d0e0f0",
                framealpha=0.8, loc="upper right")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.df); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)

    # ── Özet sekmesi ───────────────────────────────────────────────────────────
    def _show_summary(self):
        for w in self.sf.winfo_children(): w.destroy()
        params_list=[
            ("metered","Metered (mg)"),("delivered","Delivered (mg)"),
            ("fpd","FPD (mg)"),("fpf","FPF (%)"),
            ("mmad","MMAD (\u00b5m)"),("gsd","GSD"),
            ("slope","Slope"),("intercept","Intercept"),("r2","R\u00b2"),
        ]
        for sd in self.all_series:
            hf=ctk.CTkFrame(self.sf,fg_color=sd["color"],corner_radius=6,height=26)
            hf.pack(fill="x",pady=(8,2),padx=6); hf.pack_propagate(False)
            ctk.CTkLabel(hf,text=f"  {sd['name']}  |  {sd['flow']} L/min",
                font=ctk.CTkFont(size=12,weight="bold"),text_color="white",anchor="w"
                ).pack(side="left",padx=10)
            # Başlık satırı
            hdr=ctk.CTkFrame(self.sf,fg_color="#0a1422",corner_radius=4)
            hdr.pack(fill="x",padx=6,pady=1)
            ctk.CTkLabel(hdr,text="Parametre",width=180,anchor="w",
                font=ctk.CTkFont(size=10,weight="bold"),text_color="#FFC600"
                ).pack(side="left",padx=6)
            for ri in range(RUNS_PER_SERIES):
                ctk.CTkLabel(hdr,text=f"Run {ri+1}",width=85,anchor="center",
                    font=ctk.CTkFont(size=10,weight="bold"),
                    text_color=sd["color"]).pack(side="left",padx=2)
            for lbl in [self.T["mean"],self.T["sd"],self.T["rsd"]]:
                ctk.CTkLabel(hdr,text=lbl,width=80,anchor="center",
                    font=ctk.CTkFont(size=10,weight="bold"),
                    text_color="#d0e0f0").pack(side="left",padx=2)
            hlk={"fpd","fpf","mmad","gsd"}
            for j,(key,lbl) in enumerate(params_list):
                fr=ctk.CTkFrame(self.sf,
                    fg_color="#090c14" if j%2==0 else "#141824",corner_radius=3)
                fr.pack(fill="x",padx=6,pady=1)
                ctk.CTkLabel(fr,text=lbl,width=180,anchor="w",
                    font=ctk.CTkFont(size=10),text_color="#7a9abf"
                    ).pack(side="left",padx=6,pady=3)
                for run in sd["runs"]:
                    v=run.get(key,0)
                    fmt=f"{float(v):.4f}" if isinstance(v,float) else str(v)
                    ctk.CTkLabel(fr,text=fmt,width=85,anchor="center",
                        font=ctk.CTkFont(size=10),
                        text_color="#FFC600" if key in hlk else "#d0e0f0"
                        ).pack(side="left",padx=2)
                if sd["avg"] and key in sd["avg"]["params"]:
                    m_,s_,rsd_=sd["avg"]["params"][key]
                    for sv in [m_,s_,rsd_]:
                        ctk.CTkLabel(fr,text=f"{sv:.4f}",width=80,anchor="center",
                            font=ctk.CTkFont(size=10),
                            text_color="#FFC600" if key in hlk else "#f0f0f0"
                            ).pack(side="left",padx=2)
                else:
                    for _ in range(3):
                        ctk.CTkLabel(fr,text="—",width=80,anchor="center",
                            font=ctk.CTkFont(size=10),text_color="#666"
                            ).pack(side="left",padx=2)

    # ── Karşılaştırma grafiği ──────────────────────────────────────────────────
    def _plot_compare(self):
        for w in self.cf.winfo_children(): w.destroy()
        if not self.all_series: return
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        # 2x2 grid: MMAD, GSD, FPD, FPF
        params=[("mmad","MMAD (\u00b5m)"),("gsd","GSD"),
                ("fpd","FPD (mg)"),("fpf","FPF (%)")]
        for pi,(key,lbl) in enumerate(params):
            ax=fig.add_subplot(2,2,pi+1); ax.set_facecolor("#0e1525")
            names=[]; means=[]; errs=[]; colors=[]
            for sd in self.all_series:
                if sd["avg"] and key in sd["avg"]["params"]:
                    m,s,_=sd["avg"]["params"][key]
                    names.append(sd["name"]); means.append(m)
                    errs.append(s); colors.append(sd["color"])
            if not names: continue
            xs=range(len(names))
            bars=ax.bar(xs,means,color=colors,alpha=0.85,width=0.6,zorder=3)
            ax.errorbar(xs,means,yerr=errs,fmt='none',
                color='white',capsize=4,lw=1.5,zorder=4)
            for i,bar in enumerate(bars):
                ax.text(bar.get_x()+bar.get_width()/2,
                    bar.get_height()+errs[i]+max(means)*0.01,
                    f"{means[i]:.3f}",ha='center',va='bottom',
                    fontsize=7,color='white')
            ax.set_xticks(list(xs))
            ax.set_xticklabels(names,rotation=25,ha="right",fontsize=8,color="#ccc")
            ax.set_title(lbl,color="#FFC600",fontsize=10,fontweight="bold")
            ax.tick_params(colors="#7090b0",labelsize=8)
            ax.spines[:].set_color("#2a4060")
            ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
        fig.suptitle("Seri Kar\u015f\u0131la\u015ft\u0131rmas\u0131  (Ortalama ± SD)",
            color="#FFC600",fontsize=12,fontweight="bold",y=1.02)
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.cf); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)

    # ── PDF ────────────────────────────────────────────────────────────────────
    def _export_pdf(self):
        if not self.all_series:
            messagebox.showwarning("",self.T["err_nodata"]); return
        path=filedialog.asksaveasfilename(defaultextension=".pdf",
            filetypes=[("PDF","*.pdf")],
            initialfile=f"NGI_Report_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf")
        if not path: return
        meta={k:v.get() for k,v in self.mv.items()}
        try:
            from ngi_pdf import make_pdf_multi
            make_pdf_multi(path,self.all_series,meta,int(self.var_flow.get()),self.T)
            self.lbl_status.configure(text=f"PDF: {os.path.basename(path)}")
            if messagebox.askyesno("PDF","A\u00e7\u0131ls\u0131n m\u0131?"):
                if os.name=="nt": os.startfile(path)
                else: os.system(f"xdg-open '{path}'")
        except Exception as e:
            messagebox.showerror("PDF",str(e))

    # ── Temizle ────────────────────────────────────────────────────────────────
    def _clear(self):
        for sw in self.series_widgets:
            for run in sw["runs"]:
                for v in run.values(): v.set("0.000")
        self.all_series=[]
        for f in [self.rf,self.pf,self.df,self.sf,self.cf]:
            for w in f.winfo_children(): w.destroy()
        self.lbl_status.configure(text=self.T["status_ready"])

    # ── Dil değiştir ───────────────────────────────────────────────────────────
    def _toggle_lang(self):
        old_T=self.T
        self.lang="EN" if self.lang=="TR" else "TR"; self.T=L[self.lang]
        self.lbl_title.configure(text=self.T["title"])
        self.lbl_sub.configure(text=self.T["subtitle"])
        self.btn_lang.configure(text=self.T["lang_btn"])
        self.btn_add_s.configure(text=self.T["add_series"])
        self.btn_del_s.configure(text=self.T["del_series"])
        self.btn_calc.configure(text=self.T["calculate"])
        self.btn_clr.configure(text=self.T["clear"])
        self.btn_pdf.configure(text=self.T["export_pdf"])
        self.lbl_vr.configure(text=self.T["valid_range"])
        self._refresh_cutoffs()
        self.lbl_status.configure(text=self.T["status_ready"])
        # Sekme adlarını çevir
        tab_keys=["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]
        for k in tab_keys:
            try: self.tabs.rename(old_T[k], self.T[k])
            except: pass
        self.tabs.set(self.T["tab_results"])
        # Seri widget paste butonlarını çevir
        for sw in self.series_widgets:
            sw["paste_btn"].configure(text=self.T["paste_btn"])

if __name__=="__main__":
    app=NGIApp(); app.mainloop()
