"""NGI Cascade Impactor Analysis Tool v3 - Ph.Eur 2.9.18 / USP <601>"""
import subprocess, sys

def _ensure_pkg(pkg, import_name=None):
    """Paket yoksa otomatik yükle"""
    name = import_name or pkg
    try:
        __import__(name)
    except ImportError:
        print(f"{pkg} yükleniyor...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", pkg, "-q",
             "--break-system-packages"],
            creationflags=0x08000000 if sys.platform=="win32" else 0
        )

_ensure_pkg("customtkinter")
_ensure_pkg("matplotlib")
_ensure_pkg("scipy")
_ensure_pkg("Pillow", "PIL")
_ensure_pkg("reportlab")
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

# ── Icon yolu ─────────────────────────────────────────────────────────────────
def resource_path(rel):
    base = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, rel)

# ── NGI Cut-off D50 (µm) ──────────────────────────────────────────────────────
NGI_CUTOFFS = {
    15:  {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":8.61,"S3":5.39,"S4":3.30,"S5":2.08,"S6":1.36,"S7":0.98,"MOC":0.54},
    30:  {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":4.46,"S3":2.82,"S4":1.66,"S5":0.94,"S6":0.55,"S7":0.34,"MOC":0.21},
    45:  {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":3.42,"S3":2.09,"S4":1.21,"S5":0.72,"S6":0.40,"S7":0.24,"MOC":0.13},
    60:  {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":8.06,"S3":4.46,"S4":2.82,"S5":1.66,"S6":0.94,"S7":0.55,"MOC":0.34},
    90:  {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":2.08,"S3":1.36,"S4":0.98,"S5":0.55,"S6":0.34,"S7":0.21,"MOC":0.10},
    100: {"Device":999.0,"Throat":999.0,"Presep":999.0,"S1":999.0,
          "S2":1.78,"S3":1.12,"S4":0.69,"S5":0.43,"S6":0.25,"S7":0.14,"MOC":0.08},
}
STAGE_ORDER = ["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7"]
CP = ["#2E75B6","#ED7D31","#70AD47","#E84040","#7030A0","#00B0F0","#D4A000","#C00000"]

L = {
"TR":{"title":"NGI Kaskadit \u0130mpaktor Analiz Arac\u0131",
 "subtitle":"Ph.Eur 2.9.18 / USP <601>  |  Next Generation Impactor",
 "lang_btn":"English","product":"\u00dcr\u00fcn Ad\u0131","batch":"Lot No.",
 "operator":"Analist","date":"Tarih","flow_rate":"Ak\u0131\u015f H\u0131z\u0131",
 "cutoff_title":"Cut-off D50 De\u011flerleri (\u00b5m)",
 "run_label":"Run","add_run":"+ Run Ekle","del_run":"Run Sil",
 "calculate":"Hesapla","clear":"Temizle","export_pdf":"PDF Rapor",
 "tab_results":"Sonu\u00e7lar","tab_plot":"Log-Probit Graf",
 "tab_dist":"Da\u011f\u0131l\u0131m Graf","tab_summary":"\u00d6zet",
 "metered":"Metered Dose (mg)","delivered":"Delivered Dose (mg)",
 "fp_dose":"Fine Particle Dose \u22645\u00b5m (mg)","fp_frac":"Fine Particle Fraction (%)",
 "mmad":"MMAD (\u00b5m)","gsd":"GSD","r2":"R\u00b2","n_pts":"n",
 "slope":"Slope","intercept":"Intercept","mean":"Ortalama","sd":"SD","rsd":"RSD (%)",
 "valid_range":"Ge\u00e7erlilik Aral\u0131\u011f\u0131 (%)","normalize_by":"Normalize",
 "norm_ism":"ISM (Throat hari\u00e7)","norm_total":"Total (Throat dahil)",
 "settings_hdr":"Ayarlar","status_ready":"Haz\u0131r.",
 "status_done":"Hesaplama tamamland\u0131.","err_norun":"En az 1 run giriniz.",
 "err_val":"Ge\u00e7ersiz de\u011fer.","insufficient":"Yetersiz nokta (min 2)",
 "regression_hdr":"Log-Probit Regresyon",
 "x_logd":"log\u2081\u2080(D50, \u00b5m)","y_probit":"Probit z",
 "x_stage":"Stage","y_mass":"K\u00fctle (mg/at\u0131\u015f)",
 "plot_title":"Log-Probit Grafi\u011fi","dist_title":"Stage K\u00fctle Da\u011f\u0131l\u0131m\u0131",
 "pdf_saved":"PDF kaydedildi:","summary_hdr":"\u00d6zet \u2013 T\u00fcm Runlar"},
"EN":{"title":"NGI Cascade Impactor Analysis Tool",
 "subtitle":"Ph.Eur 2.9.18 / USP <601>  |  Next Generation Impactor",
 "lang_btn":"T\u00fcrk\u00e7e","product":"Product Name","batch":"Batch No.",
 "operator":"Analyst","date":"Date","flow_rate":"Flow Rate",
 "cutoff_title":"Cut-off D50 Values (\u00b5m)",
 "run_label":"Run","add_run":"+ Add Run","del_run":"Delete Run",
 "calculate":"Calculate","clear":"Clear","export_pdf":"PDF Report",
 "tab_results":"Results","tab_plot":"Log-Probit Plot",
 "tab_dist":"Distribution Plot","tab_summary":"Summary",
 "metered":"Metered Dose (mg)","delivered":"Delivered Dose (mg)",
 "fp_dose":"Fine Particle Dose \u22645\u00b5m (mg)","fp_frac":"Fine Particle Fraction (%)",
 "mmad":"MMAD (\u00b5m)","gsd":"GSD","r2":"R\u00b2","n_pts":"n",
 "slope":"Slope","intercept":"Intercept","mean":"Mean","sd":"SD","rsd":"RSD (%)",
 "valid_range":"Valid Range (%)","normalize_by":"Normalize by",
 "norm_ism":"ISM (excl. Throat)","norm_total":"Total (incl. Throat)",
 "settings_hdr":"Settings","status_ready":"Ready.",
 "status_done":"Calculation complete.","err_norun":"Enter at least 1 run.",
 "err_val":"Invalid value.","insufficient":"Insufficient points (min 2)",
 "regression_hdr":"Log-Probit Regression",
 "x_logd":"log\u2081\u2080(D50, \u00b5m)","y_probit":"Probit z",
 "x_stage":"Stage","y_mass":"Mass (mg/act.)",
 "plot_title":"Log-Probit Plot","dist_title":"Stage Mass Distribution",
 "pdf_saved":"PDF saved:","summary_hdr":"Summary \u2013 All Runs"},
}

# ── Hesaplama ─────────────────────────────────────────────────────────────────
def calc_run(masses, flow, lo=15, hi=85, use_ism=True):
    co = NGI_CUTOFFS[flow]
    device = masses.get("Device", 0)
    throat = masses.get("Throat", 0)
    presep = masses.get("Presep", 0)
    metered  = sum(masses.values())
    # ISM = S1..MOC (Device, Throat, Presep hariç)
    ism = sum(masses.get(s,0) for s in ["S1","S2","S3","S4","S5","S6","S7","MOC"])
    delivered = ism  # Emitted = ISM
    denom = ism if use_ism else metered
    ism_stages = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
    cum = []
    for s in STAGE_ORDER:
        mass = masses.get(s, 0)
        u = 0.0
        if s in ism_stages and denom > 0:
            # CITDAS: inclusive undersize (o stage dahil)
            u = sum(masses.get(x,0) for x in ism_stages
                    if co.get(x,0) <= co.get(s,0)) / denom * 100
        cum.append({"stage":s,"d50":co.get(s,0),"mass":mass,"u_pct":u})
    valid = [r for r in cum if r["stage"]!="Presep" and lo<r["u_pct"]<hi]
    res = {"metered":metered,"delivered":delivered,"cum_data":cum,
           "valid":valid,"masses":masses,"flow":flow}
    if len(valid)<2:
        res["error"]="insufficient"; res["n"]=len(valid); return res
    x = np.array([math.log10(v["d50"]) for v in valid])
    y = np.array([norm.ppf(v["u_pct"]/100) for v in valid])
    b = np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
    a = y.mean()-b*x.mean()
    yp = a+b*x; ss_r=np.sum((y-yp)**2); ss_t=np.sum((y-y.mean())**2)
    r2 = 1-ss_r/ss_t if ss_t>0 else 1.0
    z5 = a+b*math.log10(5.0)
    res.update({"n":len(valid),"a":a,"b":b,"slope":b,"intercept":a+5,"r2":r2,
                "mmad":10**(-a/b),"gsd":10**(1/b),
                "fpd":norm.cdf(z5)*ism,"fpf":norm.cdf(z5)*ism/metered*100,"x_reg":x,"y_reg":y})
    return res

def calc_summary(runs):
    out={}
    for p in ["mmad","gsd","fpd","fpf","metered","delivered","slope","intercept","r2","n"]:
        vals=[r[p] for r in runs if "error" not in r and p in r]
        if not vals: out[p]=(None,None,None); continue
        m=float(np.mean(vals)); s=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
        out[p]=(m, s, s/m*100 if m else 0.0)
    return out

# ── PDF (Türkçe DejaVu, siyah-beyaz, KeepTogether) ───────────────────────────
def make_pdf(path, run_results, meta, flow, T):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                    Paragraph, Spacer, Image,
                                    HRFlowable)
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import io as _io

    # Türkçe font kaydı
    # Uygulama dizinindeki font
    _app_dir = os.path.dirname(os.path.abspath(__file__))
    font_paths = [
        os.path.join(_app_dir, 'DejaVuSans.ttf'),
        '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
        'C:/Windows/Fonts/DejaVuSans.ttf',
        'C:/Windows/Fonts/arial.ttf',
        'C:/Windows/Fonts/calibri.ttf',
    ]
    font_bold_paths = [
        os.path.join(_app_dir, 'DejaVuSans-Bold.ttf'),
        '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
        'C:/Windows/Fonts/DejaVuSans-Bold.ttf',
        'C:/Windows/Fonts/arialbd.ttf',
        'C:/Windows/Fonts/calibrib.ttf',
    ]
    fn = "DVSans"; fnb = "DVSansBold"
    registered = False
    for fp, fb in zip(font_paths, font_bold_paths):
        if os.path.exists(fp) and os.path.exists(fb):
            try:
                pdfmetrics.registerFont(TTFont(fn, fp))
                pdfmetrics.registerFont(TTFont(fnb, fb))
                registered = True
                break
            except: pass
    if not registered:
        fn = "Helvetica"; fnb = "Helvetica-Bold"

    # Siyah-beyaz renk paleti
    BLK  = colors.black
    DGR  = colors.HexColor("#404040")   # Koyu gri (başlık)
    MGR  = colors.HexColor("#555555")   # Orta gri
    LGR  = colors.HexColor("#AAAAAA")   # Açık gri
    XLGR = colors.HexColor("#E8E8E8")   # Çok açık gri (satır bg)
    XXLGR= colors.HexColor("#F5F5F5")
    WHT  = colors.white
    HLW  = colors.HexColor("#EEEEEE")   # Highlight (sarı yerine açık gri)

    W_ = A4[0]-3.6*cm

    def ps(name, sz=8, bold=False, clr=BLK, align=TA_LEFT, leading=None):
        return ParagraphStyle(name, fontSize=sz, leading=leading or sz*1.1,
            fontName=fnb if bold else fn, textColor=clr, alignment=align,
            splitLongWords=False, wordWrap='LTR',
            spaceBefore=0, spaceAfter=0)

    sN  = ps("n")
    sB  = ps("b",sz=7.5,bold=True)
    sT  = ps("t",sz=15,bold=True,clr=BLK,align=TA_CENTER)
    sSb = ps("sb",sz=9,clr=MGR,align=TA_CENTER)
    sH  = ps("h",sz=9,bold=True,clr=BLK)
    sL  = ps("l",sz=7,clr=MGR)
    sV  = ps("v",sz=9,bold=True)
    sF  = ps("f",sz=7,clr=MGR,align=TA_CENTER)

    def sec_tbl(text, style, bg):
        t = Table([[Paragraph(text, style)]], colWidths=[W_])
        t.setStyle(TableStyle([
            ("TOPPADDING",(0,0),(-1,-1),1),("BOTTOMPADDING",(0,0),(-1,-1),1),
            ("LEFTPADDING",(0,0),(-1,-1),0),("RIGHTPADDING",(0,0),(-1,-1),0),
            ("LINEBELOW",(0,0),(-1,-1),0.8,MGR),
        ]))
        return t

    def grid_style(hdr_bg=DGR, alt_bg=XLGR, hl_cols=None):
        ts = [
            ("TEXTCOLOR",(0,0),(-1,-1),BLK),
            ("FONTNAME",(0,0),(-1,0),fnb),
            ("FONTSIZE",(0,0),(-1,-1),8),
            ("TOPPADDING",(0,0),(-1,-1),1),
            ("BOTTOMPADDING",(0,0),(-1,-1),1),
            ("LEFTPADDING",(0,0),(-1,-1),0),
            ("RIGHTPADDING",(0,0),(-1,-1),0),
            ("LEFTPADDING",(0,0),(-1,-1),0),
            ("RIGHTPADDING",(0,0),(-1,-1),0),
            ("ALIGN",(0,0),(-1,-1),"CENTER"),
            ("BOX",(0,0),(-1,-1),0.8,MGR),
            ("INNERGRID",(0,0),(-1,-1),0.3,LGR),
            ("LINEBELOW",(0,0),(-1,0),1.0,MGR),
        ]
        return ts

    doc = SimpleDocTemplate(path, pagesize=A4,
        topMargin=1.5*cm, bottomMargin=1.5*cm,
        leftMargin=1.8*cm, rightMargin=1.8*cm)
    story = []

    # Başlık
    story.append(sec_tbl(T["title"], sT, DGR))
    story.append(sec_tbl(T["subtitle"], sSb, MGR))
    story.append(Spacer(1,0.25*cm))

    # Meta
    co = NGI_CUTOFFS[flow]
    mr = [[Paragraph(T["product"],sL),Paragraph(meta.get("product",""),sV),
           Paragraph(T["batch"],sL),Paragraph(meta.get("batch",""),sV),
           Paragraph(T["flow_rate"],sL),Paragraph(f"{flow} L/min",sV)],
          [Paragraph(T["operator"],sL),Paragraph(meta.get("operator",""),sV),
           Paragraph(T["date"],sL),Paragraph(meta.get("date",""),sV),
           Paragraph("",sL),Paragraph("",sV)]]
    mt = Table(mr, colWidths=[W_/6]*6)
    mt.setStyle(TableStyle([("BOX",(0,0),(-1,-1),0.5,MGR),("INNERGRID",(0,0),(-1,-1),0.3,LGR),
        ("TOPPADDING",(0,0),(-1,-1),1),("BOTTOMPADDING",(0,0),(-1,-1),1),
        ("LEFTPADDING",(0,0),(-1,-1),0),("RIGHTPADDING",(0,0),(-1,-1),0)]))
    story.append(mt); story.append(Spacer(1,0.2*cm))

    # Cut-off tablosu
    co_h = [Paragraph("D50\n(\u00b5m)",sB)] + [Paragraph(s,sB) for s in STAGE_ORDER]
    co_v = [Paragraph(f"{flow} L/min",sN)] + [Paragraph(f"{co[s]:.2f}",sN) for s in STAGE_ORDER]
    ct = Table([co_h,co_v], colWidths=[W_/9]*9)
    ts = grid_style()

    ct.setStyle(TableStyle(ts))
    story.append(ct); story.append(Spacer(1,0.3*cm))

    valid_runs = [r for r in run_results if "error" not in r]

    for res in run_results:
        rn = res["run_no"]
        run_elements = []

        story.append(sec_tbl(f"  {T['run_label']} {rn}", sH, DGR))
        story.append(Spacer(1,0.1*cm))

        if "error" in res:
            story.append(Paragraph(f"\u26a0  {T['insufficient']}", sB))
            story.append(Spacer(1,0.2*cm))
            continue

        # Kümülatif tablo
        story.append(sec_tbl("Measured Drug per Discharge  [mg/actuation]", sH, MGR))
        valid_s = {v["stage"] for v in res["valid"]}
        ch = [Paragraph(t,sB) for t in ["Stage","D50\n(\u00b5m)","Mass\n(mg)",
              "Cum.Mass (mg)","Cum. (%)","Valid","Probit z"]]
        crows = [ch]; cm_ = 0
        for row in res["cum_data"]:
            cm_ += row["mass"]; iv = row["stage"] in valid_s
            try: pz = f"{norm.ppf(row['u_pct']/100):.4f}" if 0<row["u_pct"]<100 else "\u2014"
            except: pz = "\u2014"
            crows.append([Paragraph(row["stage"], sB if iv else sN),
                Paragraph(f"{row['d50']:.3f}",sN), Paragraph(f"{row['mass']:.4f}",sN),
                Paragraph(f"{cm_:.4f}",sN), Paragraph(f"{row['u_pct']:.3f}",sN),
                Paragraph("Yes" if iv else "No", sN), Paragraph(pz,sN)])
        ccw = [W_*x for x in [0.1017,0.1557,0.1723,0.2450,0.1524,0.0411,0.1318]]
        cmt = Table(crows, colWidths=ccw)
        cts = grid_style()
        for i,row in enumerate(res["cum_data"],1):
            if row["stage"] in valid_s: cts.append(("FONTNAME",(0,i),(-1,i),fnb))
        cmt.setStyle(TableStyle(cts))
        story.append(cmt); story.append(Spacer(1,0.15*cm))

        # Dozaj tablosu
        story.append(sec_tbl("Particle Size Characterisation", sH, MGR))
        dh = [Paragraph(t,sB) for t in ["Metered (mg)","Delivered (mg)",
              "FPD (mg)","FPF (%)","MMAD (\u00b5m)","GSD","R\u00b2","n","Intercept","Slope"]]
        dv = [Paragraph(f"{res['metered']:.4f}",sN), Paragraph(f"{res['delivered']:.4f}",sN),
              Paragraph(f"{res['fpd']:.4f}",sB), Paragraph(f"{res['fpf']:.3f}",sB),
              Paragraph(f"{res['mmad']:.4f}",sB), Paragraph(f"{res['gsd']:.4f}",sB),
              Paragraph(f"{res['r2']:.4f}",sN), Paragraph(str(res["n"]),sN),
              Paragraph(f"{res['intercept']:.4f}",sN), Paragraph(f"{res['slope']:.4f}",sN)]
        # Kolon genişlikleri: tek satır için min genişliğe göre hesaplanmış
        dt_cw = [W_*x for x in [0.1506,0.1630,0.1068,0.0908,0.1311,0.0754,0.0754,0.0259,0.1056,0.0754]]
        dt = Table([dh,dv], colWidths=dt_cw)
        dts = grid_style()
        dts.append(("FONTNAME",(2,1),(5,1),fnb))
        dts.append(("LEFTPADDING",(7,0),(7,-1),0))
        dts.append(("RIGHTPADDING",(7,0),(7,-1),0))
        dt.setStyle(TableStyle(dts))
        story.append(dt); story.append(Spacer(1,0.25*cm))


    # Özet tablosu
    if len(valid_runs) > 1:
        from reportlab.platypus import KeepTogether as KT
        _sum_hdr = sec_tbl("Summary Statistics \u2013 All Runs", sH, DGR)
        sm = calc_summary(valid_runs)
        psum = [("metered","Metered Dose (mg)"),("delivered","Delivered Dose (mg)"),
                ("fpd","Fine Particle Dose (mg)"),("fpf","Fine Particle Fraction (%)"),
                ("mmad","MMAD (\u00b5m)"),("gsd","GSD"),("slope","Slope"),
                ("intercept","Intercept"),("r2","R\u00b2"),("n","n")]
        hlk = {"fpd","fpf","mmad","gsd"}
        sh = [Paragraph("Parameter",sB)] + \
             [Paragraph(f"Run {r['run_no']}",sB) for r in valid_runs] + \
             [Paragraph(t,sB) for t in ["Mean","SD","RSD (%)"]]
        srows = [sh]
        for key,lbl in psum:
            row = [Paragraph(lbl, sB if key in hlk else sN)]
            for r in valid_runs:
                v = r.get(key,0)
                row.append(Paragraph(f"{float(v):.4f}" if isinstance(v,float) else str(v),
                                     sB if key in hlk else sN))
            m_,s_,rsd_ = sm.get(key,(None,None,None))
            for sv in [m_,s_,rsd_]:
                row.append(Paragraph(f"{sv:.4f}" if sv is not None else "\u2014", sN))
            srows.append(row)
        nc = 1+len(valid_runs)+3
        rem = W_ - W_*0.22
        scw = [W_*0.22] + [rem/(nc-1)]*(nc-1)
        stbl = Table(srows, colWidths=scw[:nc])
        sts = grid_style()
        for i,(key,_) in enumerate(psum,1):
            if key in hlk: sts.append(("FONTNAME",(0,i),(-1,i),fnb))
        stbl.setStyle(TableStyle(sts))
        _sum_elements = [_sum_hdr, stbl, Spacer(1,0.1*cm)]
        

    # Grafikler
    import matplotlib.pyplot as plt
    import io as _io2

    fig_lp = Figure(figsize=(6.5,3.8))
    ax_lp  = fig_lp.add_subplot(111)
    ax_lp.set_facecolor("white"); fig_lp.patch.set_facecolor("white")
    markers = ['o','s','^','D','v','<','>','p']
    greys   = ["#000000","#444444","#888888","#222222","#666666","#aaaaaa"]
    for i,res in enumerate(valid_runs):
        gc = greys[i%len(greys)]
        ax_lp.scatter(res["x_reg"],res["y_reg"],color=gc,s=60,zorder=5,
                      marker=markers[i%len(markers)],label=f"Run {res['run_no']}")
        xl = np.linspace(min(res["x_reg"])-0.15,max(res["x_reg"])+0.15,100)
        ax_lp.plot(xl,res["a"]+res["b"]*xl,color=gc,lw=1.4,ls=['-','--',':','-.'][i%4])
        mx = math.log10(res["mmad"])
        ax_lp.axvline(mx,color=gc,ls="--",lw=0.8,alpha=0.7)
    ax_lp.set_xlabel("log\u2081\u2080(D50, \u00b5m)",fontsize=8)
    ax_lp.set_ylabel("Probit z",fontsize=8)
    ax_lp.set_title("Log-Probit Graph",fontsize=9,fontweight="bold")
    ax_lp.grid(True,color="#cccccc",ls="--",alpha=0.6)
    ax_lp.legend(fontsize=7); ax_lp.tick_params(labelsize=7)
    fig_lp.tight_layout()
    buf_lp = _io2.BytesIO()
    fig_lp.savefig(buf_lp,format="png",dpi=150,bbox_inches="tight",facecolor="white")
    buf_lp.seek(0); plt.close('all')

    stages = STAGE_ORDER+["MOC"]
    fig_d = Figure(figsize=(5.5,3.8))
    ax_d  = fig_d.add_subplot(111)
    ax_d.set_facecolor("white"); fig_d.patch.set_facecolor("white")
    nv=len(valid_runs); bw=0.8/max(nv,1); xp=np.arange(len(stages))
    hatches = ['','//','xx','..','\\\\','||','--','++']
    for i,res in enumerate(valid_runs):
        gc=greys[i%len(greys)]; ms=[res["masses"].get(s,0) for s in stages]
        off=(i-nv/2+0.5)*bw
        ax_d.bar(xp+off,ms,width=bw*0.88,color=gc,alpha=0.85,
                 hatch=hatches[i%len(hatches)],edgecolor="white",
                 label=f"Run {res['run_no']}")
    ax_d.set_xticks(xp); ax_d.set_xticklabels(stages,rotation=30,ha="right",fontsize=7)
    ax_d.set_xlabel("Stage",fontsize=8); ax_d.set_ylabel("Mass (mg)",fontsize=8)
    ax_d.set_title(f"Drug Distribution  [{flow} L/min]",fontsize=9,fontweight="bold")
    ax_d.grid(True,axis="y",color="#cccccc",ls="--",alpha=0.6)
    ax_d.legend(fontsize=7); ax_d.tick_params(labelsize=7)
    fig_d.tight_layout()
    buf_d = _io2.BytesIO()
    fig_d.savefig(buf_d,format="png",dpi=150,bbox_inches="tight",facecolor="white")
    buf_d.seek(0); plt.close('all')

    img_lp = Image(buf_lp, width=W_*0.54, height=W_*0.30)
    img_d  = Image(buf_d,  width=W_*0.44, height=W_*0.30)
    gt = Table([[img_lp,img_d]], colWidths=[W_*0.56,W_*0.44])
    gt.setStyle(TableStyle([("ALIGN",(0,0),(-1,-1),"CENTER"),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),("BOX",(0,0),(-1,-1),0.5,LGR)]))
    # Summary + Grafikler aynı sayfada
    if len(valid_runs) > 1:
        _sum_elements.append(gt)
        _sum_elements.append(Spacer(1,0.2*cm))
        story.append(KT(_sum_elements))
    else:
        story.append(gt)
        story.append(Spacer(1,0.2*cm))
    story.append(HRFlowable(width=W_,thickness=0.5,color=MGR))
    story.append(Paragraph(
        f"NGI Analysis Tool  |  Ph.Eur 2.9.18 / USP <601>  |  {datetime.now().strftime('%d.%m.%Y %H:%M')}",
        sF))
    doc.build(story)


# ── Ana Uygulama ──────────────────────────────────────────────────────────────
class NGIApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.lang="TR"; self.T=L["TR"]; self.run_results=[]
        ctk.set_appearance_mode("dark"); ctk.set_default_color_theme("blue")
        self.title(self.T["title"]); self.geometry("1380x900"); self.minsize(1100,720)
        # İkon
        ico = resource_path("icon.ico")
        if os.path.exists(ico):
            try: self.iconbitmap(ico)
            except: pass
        self._ui()
        self.lbl_status.configure(text=self.T["status_ready"])

    def _ui(self):
        hdr=ctk.CTkFrame(self,fg_color="#002D62",corner_radius=0,height=56)
        hdr.pack(fill="x"); hdr.pack_propagate(False)
        self.lbl_title=ctk.CTkLabel(hdr,text=self.T["title"],
            font=ctk.CTkFont(size=16,weight="bold"),text_color="#FFC600")
        self.lbl_title.pack(side="left",padx=16,pady=6)
        self.lbl_sub=ctk.CTkLabel(hdr,text=self.T["subtitle"],
            font=ctk.CTkFont(size=10),text_color="#aac8e8")
        self.lbl_sub.pack(side="left",padx=4)
        self.btn_lang=ctk.CTkButton(hdr,text=self.T["lang_btn"],width=85,height=30,
            command=self._lang,fg_color="#001a40",hover_color="#003580")
        self.btn_lang.pack(side="right",padx=14)
        body=ctk.CTkFrame(self,fg_color="transparent"); body.pack(fill="both",expand=True)
        self.left=ctk.CTkScrollableFrame(body,width=430,fg_color="#141824",corner_radius=0)
        self.left.pack(side="left",fill="y")
        self._left()
        right=ctk.CTkFrame(body,fg_color="#0e1219"); right.pack(side="left",fill="both",expand=True)
        self._right(right)
        sb=ctk.CTkFrame(self,height=26,fg_color="#090c12",corner_radius=0)
        sb.pack(fill="x",side="bottom"); sb.pack_propagate(False)
        self.lbl_status=ctk.CTkLabel(sb,text="",anchor="w",
            font=ctk.CTkFont(size=11),text_color="#4a7a9a")
        self.lbl_status.pack(side="left",padx=12)

    def _sec(self,t):
        ctk.CTkLabel(self.left,text=t,anchor="w",
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFC600").pack(fill="x",padx=8,pady=(10,2))

    def _card(self):
        f=ctk.CTkFrame(self.left,fg_color="#1c2336",corner_radius=8); f.pack(fill="x",padx=6,pady=3)
        i=ctk.CTkFrame(f,fg_color="transparent"); i.pack(fill="x",padx=10,pady=8); return i

    def _left(self):
        self._sec("● "+self.T["tab_results"])
        mc=self._card(); self.mv={}
        for k,d in [("product",""),("batch",""),("operator",""),
                    ("date",datetime.now().strftime("%d.%m.%Y"))]:
            r=ctk.CTkFrame(mc,fg_color="transparent"); r.pack(fill="x",pady=2)
            ctk.CTkLabel(r,text=self.T[k],width=88,anchor="w",
                font=ctk.CTkFont(size=11),text_color="#7a9abf").pack(side="left")
            v=ctk.StringVar(value=d); self.mv[k]=v
            ctk.CTkEntry(r,textvariable=v,height=26,
                font=ctk.CTkFont(size=11)).pack(side="left",fill="x",expand=True,padx=4)

        self._sec("● "+self.T["flow_rate"])
        fc=self._card()
        fr=ctk.CTkFrame(fc,fg_color="transparent"); fr.pack(fill="x",pady=2)
        ctk.CTkLabel(fr,text=self.T["flow_rate"],width=100,anchor="w",
            font=ctk.CTkFont(size=11),text_color="#7a9abf").pack(side="left")
        self.var_flow=ctk.StringVar(value="30")
        ctk.CTkOptionMenu(fr,values=["15","30","45","60","90","100"],
            variable=self.var_flow,command=self._on_flow,width=90,height=28,
            font=ctk.CTkFont(size=12,weight="bold")).pack(side="left",padx=6)
        ctk.CTkLabel(fr,text="L/min",font=ctk.CTkFont(size=11),
            text_color="#7a9abf").pack(side="left")
        self.cbox=ctk.CTkFrame(fc,fg_color="#090c14",corner_radius=6)
        self.cbox.pack(fill="x",pady=(6,2)); self._co()

        self._sec("● "+self.T["settings_hdr"])
        sc=self._card()
        r1=ctk.CTkFrame(sc,fg_color="transparent"); r1.pack(fill="x",pady=2)
        self.lbl_vr=ctk.CTkLabel(r1,text=self.T["valid_range"],width=150,anchor="w",
            font=ctk.CTkFont(size=11),text_color="#7a9abf"); self.lbl_vr.pack(side="left")
        self.var_lo=ctk.StringVar(value="15"); self.var_hi=ctk.StringVar(value="85")
        ctk.CTkEntry(r1,textvariable=self.var_lo,width=50,height=26).pack(side="left",padx=2)
        ctk.CTkLabel(r1,text="\u2013").pack(side="left",padx=2)
        ctk.CTkEntry(r1,textvariable=self.var_hi,width=50,height=26).pack(side="left",padx=2)
        r2=ctk.CTkFrame(sc,fg_color="transparent"); r2.pack(fill="x",pady=2)
        self.lbl_nb=ctk.CTkLabel(r2,text=self.T["normalize_by"],width=150,anchor="w",
            font=ctk.CTkFont(size=11),text_color="#7a9abf"); self.lbl_nb.pack(side="left")
        self.var_norm=ctk.StringVar(value="ism")
        self.rb1=ctk.CTkRadioButton(r2,text=self.T["norm_ism"],variable=self.var_norm,
            value="ism",font=ctk.CTkFont(size=10)); self.rb1.pack(side="left",padx=4)
        self.rb2=ctk.CTkRadioButton(r2,text=self.T["norm_total"],variable=self.var_norm,
            value="total",font=ctk.CTkFont(size=10)); self.rb2.pack(side="left",padx=4)

        self._sec("● Stage Masses (mg/actuation)")
        self.rbox=ctk.CTkFrame(self.left,fg_color="transparent")
        self.rbox.pack(fill="x",padx=4,pady=2)
        self.rw=[]; self._add()

        br=ctk.CTkFrame(self.left,fg_color="transparent"); br.pack(fill="x",padx=6,pady=4)
        self.btn_add=ctk.CTkButton(br,text=self.T["add_run"],command=self._add,
            height=30,width=115,fg_color="#1a4d30",hover_color="#256640")
        self.btn_add.pack(side="left",padx=3)
        self.btn_del=ctk.CTkButton(br,text=self.T["del_run"],command=self._del,
            height=30,width=100,fg_color="#4d1a1a",hover_color="#662525")
        self.btn_del.pack(side="left",padx=3)

        ar=ctk.CTkFrame(self.left,fg_color="transparent"); ar.pack(fill="x",padx=6,pady=8)
        self.btn_calc=ctk.CTkButton(ar,text=self.T["calculate"],command=self._calc,
            height=38,width=140,font=ctk.CTkFont(size=13,weight="bold"),
            fg_color="#002D62",hover_color="#003a80")
        self.btn_calc.pack(side="left",padx=3)
        self.btn_clear=ctk.CTkButton(ar,text=self.T["clear"],command=self._clr,
            height=38,width=90,fg_color="#333",hover_color="#555")
        self.btn_clear.pack(side="left",padx=3)
        self.btn_pdf=ctk.CTkButton(ar,text=self.T["export_pdf"],command=self._pdf,
            height=38,width=120,fg_color="#1a4d30",hover_color="#256640")
        self.btn_pdf.pack(side="left",padx=3)

    def _right(self,parent):
        self.tabs=ctk.CTkTabview(parent,fg_color="#0e1219",
            segmented_button_fg_color="#141824",
            segmented_button_selected_color="#002D62",
            segmented_button_unselected_color="#1c2336")
        self.tabs.pack(fill="both",expand=True,padx=4,pady=4)
        for k in ["tab_results","tab_plot","tab_dist","tab_summary"]:
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

    def _co(self):
        for w in self.cbox.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        ctk.CTkLabel(self.cbox,
            text=f" {self.T['cutoff_title']}  [{flow} L/min]",
            font=ctk.CTkFont(size=10,weight="bold"),
            text_color="#FFC600",anchor="w").pack(fill="x",padx=6,pady=(4,2))
        gf=ctk.CTkFrame(self.cbox,fg_color="transparent"); gf.pack(padx=6,pady=(0,6))
        for s in [s for s in STAGE_ORDER if s not in ("Device","Throat","Presep","S1")] + ["MOC"]:
            c=ctk.CTkFrame(gf,fg_color="#001a40",corner_radius=4,width=48)
            c.pack(side="left",padx=2); c.pack_propagate(False)
            ctk.CTkLabel(c,text=s.replace("Presep","Pre"),
                font=ctk.CTkFont(size=8),text_color="#aac8e8").pack(pady=(2,0))
            ctk.CTkLabel(c,text=f"{co[s]:.2f}",
                font=ctk.CTkFont(size=10,weight="bold"),
                text_color="#FFC600").pack(pady=(0,2))

    def _on_flow(self,v): self._co()

    def _add(self):
        idx=len(self.rw)+1
        frame=ctk.CTkFrame(self.rbox,fg_color="#1c2336",corner_radius=8)
        frame.pack(fill="x",pady=3,padx=2)
        hf=ctk.CTkFrame(frame,fg_color="#001a40",corner_radius=6)
        hf.pack(fill="x",padx=6,pady=(6,2))
        ctk.CTkLabel(hf,text=f"  {self.T['run_label']} {idx}",
            font=ctk.CTkFont(size=11,weight="bold"),
            text_color=CP[(idx-1)%len(CP)]).pack(side="left",pady=4)
        entries={}
        # Satır 1: Device, Throat, Presep  Satır 2: S1-S4  Satır 3: S5-S7+MOC
        for row_s in [["Device","Throat","Presep"],
                      STAGE_ORDER[3:7],
                      STAGE_ORDER[7:]+["MOC"]]:
            rf=ctk.CTkFrame(frame,fg_color="transparent"); rf.pack(fill="x",padx=6,pady=2)
            for s in row_s:
                cf=ctk.CTkFrame(rf,fg_color="transparent")
                cf.pack(side="left",expand=True,fill="x",padx=2)
                ctk.CTkLabel(cf,text=s,
                    font=ctk.CTkFont(size=9),text_color="#5a8ab0",anchor="center").pack()
                v=ctk.StringVar(value="0.000")
                ctk.CTkEntry(cf,textvariable=v,height=26,width=58,
                    font=ctk.CTkFont(size=11),justify="center").pack()
                entries[s]=v
        self.rw.append({"frame":frame,"entries":entries})

    def _del(self):
        if len(self.rw)>1: self.rw.pop()["frame"].destroy()

    def _calc(self):
        try: lo=float(self.var_lo.get()); hi=float(self.var_hi.get())
        except: lo,hi=15,85
        flow=int(self.var_flow.get()); use_ism=self.var_norm.get()=="ism"
        self.run_results=[]
        for i,rw in enumerate(self.rw):
            try:
                all_keys = list(STAGE_ORDER) + ["MOC"]
                m={s:float(rw["entries"][s].get().replace(",",".")) for s in all_keys if s in rw["entries"]}
            except:
                messagebox.showerror("Error",f"{self.T['err_val']} (Run {i+1})"); return
            res=calc_run(m,flow,lo,hi,use_ism)
            res["run_no"]=i+1; self.run_results.append(res)
        self._res()
        self.after(150, self._lp)
        self.after(300, self._dist)
        self.after(450, self._sum)
        self.tabs.set(self.T["tab_results"])
        self.lbl_status.configure(text=self.T["status_done"])

    def _res(self):
        for w in self.rf.winfo_children(): w.destroy()
        flow=int(self.var_flow.get())
        for res in self.run_results:
            rn=res["run_no"]; col=CP[(rn-1)%len(CP)]
            hf=ctk.CTkFrame(self.rf,fg_color=col,corner_radius=6,height=30)
            hf.pack(fill="x",pady=(8,2),padx=6); hf.pack_propagate(False)
            ctk.CTkLabel(hf,text=f"  {self.T['run_label']} {rn}  \u2014  {flow} L/min",
                font=ctk.CTkFont(size=12,weight="bold"),text_color="white",
                anchor="w").pack(side="left",padx=10)
            if "error" in res:
                ctk.CTkLabel(self.rf,text=f"  \u26a0  {self.T['insufficient']}",
                    text_color="#ff6644",anchor="w").pack(fill="x",padx=12,pady=4)
                continue
            self._rc(self.rf,"📊  "+self.T["tab_results"],
                [(self.T["metered"],f"{res['metered']:.4f} mg"),
                 (self.T["delivered"],f"{res['delivered']:.4f} mg"),
                 (self.T["fp_dose"],f"{res['fpd']:.4f} mg",True),
                 (self.T["fp_frac"],f"{res['fpf']:.3f} %",True)])
            self._rc(self.rf,"📐  MMAD / GSD",
                [(self.T["mmad"],f"{res['mmad']:.4f} \u00b5m",True),
                 (self.T["gsd"],f"{res['gsd']:.4f}",True),
                 (self.T["r2"],f"{res['r2']:.4f}"),
                 (self.T["n_pts"],str(res["n"]))])
            self._rc(self.rf,"📈  "+self.T["regression_hdr"],
                [(self.T["slope"],f"{res['slope']:.4f}"),
                 (self.T["intercept"],f"{res['intercept']:.4f}"),
                 (self.T["r2"],f"{res['r2']:.4f}"),
                 (self.T["n_pts"],str(res["n"]))])
            self._ct(res)

    def _rc(self,parent,title,rows):
        f=ctk.CTkFrame(parent,fg_color="#1c2336",corner_radius=8)
        f.pack(fill="x",padx=6,pady=3)
        ctk.CTkLabel(f,text=title,font=ctk.CTkFont(size=11,weight="bold"),
            text_color="#FFC600",anchor="w").pack(fill="x",padx=10,pady=(6,2))
        inn=ctk.CTkFrame(f,fg_color="transparent"); inn.pack(fill="x",padx=8,pady=(0,6))
        for i,row in enumerate(rows):
            hl=len(row)>2 and row[2]
            fr=ctk.CTkFrame(inn,fg_color="#090c14" if i%2==0 else "transparent",corner_radius=3)
            fr.pack(fill="x",pady=1)
            ctk.CTkLabel(fr,text=row[0],width=230,anchor="w",
                font=ctk.CTkFont(size=11),text_color="#7a9abf").pack(side="left",padx=8,pady=3)
            ctk.CTkLabel(fr,text=row[1],anchor="e",
                font=ctk.CTkFont(size=11,weight="bold" if hl else "normal"),
                text_color="#FFC600" if hl else "#e0e0e0").pack(side="right",padx=12,pady=3)

    def _ct(self,res):
        f=ctk.CTkFrame(self.rf,fg_color="#1c2336",corner_radius=8)
        f.pack(fill="x",padx=6,pady=3)
        ctk.CTkLabel(f,text="📋  Cumulative Table",
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFC600",
            anchor="w").pack(fill="x",padx=10,pady=(6,2))
        inn=ctk.CTkFrame(f,fg_color="transparent"); inn.pack(fill="x",padx=8,pady=(0,6))
        hdr=ctk.CTkFrame(inn,fg_color="#001a40",corner_radius=4); hdr.pack(fill="x",pady=1)
        for txt,w in [("Stage",80),("D50",70),("Mass",80),("Cum %",115),("Valid",48),("Probit z",85)]:
            ctk.CTkLabel(hdr,text=txt,width=w,anchor="center",
                font=ctk.CTkFont(size=10,weight="bold"),
                text_color="#FFC600").pack(side="left",padx=2,pady=3)
        vs={v["stage"] for v in res["valid"]}
        for i,row in enumerate(res["cum_data"]):
            iv=row["stage"] in vs
            fr=ctk.CTkFrame(inn,fg_color="#1a3050" if iv else "#090c14" if i%2==0 else "transparent",
                corner_radius=3); fr.pack(fill="x",pady=1)
            try: pz=f"{norm.ppf(row['u_pct']/100):.4f}" if 0<row["u_pct"]<100 else "\u2014"
            except: pz="\u2014"
            for txt,w,tc in [
                (row["stage"],80,"#90c8e0" if iv else "#8aabcc"),
                (f"{row['d50']:.3f}",70,"#e0e0e0"),
                (f"{row['mass']:.4f}",80,"#e0e0e0"),
                (f"{row['u_pct']:.3f}%",115,"#FFC600" if iv else "#c0d8f0"),
                ("\u2713" if iv else "\u2013",48,"#50e080" if iv else "#666"),
                (pz,85,"#c8e8ff" if iv else "#8aabcc")]:
                ctk.CTkLabel(fr,text=txt,width=w,anchor="center",
                    font=ctk.CTkFont(size=10),text_color=tc).pack(side="left",padx=2,pady=2)

    def _lp(self):
        for w in self.pf.winfo_children(): w.destroy()
        vr=[r for r in self.run_results if "error" not in r]
        if not vr: return
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        markers=['o','s','^','D','v']
        for res in vr:
            c=CP[(res["run_no"]-1)%len(CP)]; rn=res["run_no"]
            ax.scatter(res["x_reg"],res["y_reg"],color=c,s=70,zorder=5,
                marker=markers[(rn-1)%len(markers)],label=f"Run {rn}")
            xl=np.linspace(min(res["x_reg"])-0.15,max(res["x_reg"])+0.15,120)
            ax.plot(xl,res["a"]+res["b"]*xl,color=c,lw=1.6)
            mx=math.log10(res["mmad"])
            ax.axvline(mx,color=c,ls="--",lw=1,alpha=0.6)
        # MMAD etiketleri - çakışma önleme
        mmad_xs = [math.log10(r["mmad"]) for r in vr]
        used_y = []
        ymin = ax.get_ylim()[0] if ax.get_ylim()[0] > -4 else -3.0
        for i,res in enumerate(vr):
            c=CP[(res["run_no"]-1)%len(CP)]
            mx=math.log10(res["mmad"])
            y_pos = ymin + 0.1
            # Çakışma kontrolü: diğer etiketlerle mesafe
            for uy,ux in used_y:
                if abs(mx-ux) < 0.12:
                    y_pos = uy + 0.45
            used_y.append((y_pos, mx))
            ax.text(mx+0.015, y_pos,
                f"MMAD={res['mmad']:.3f}\u00b5m",
                color=c,fontsize=7.5,rotation=90,va="bottom",ha="left")

        # Sağ alt köşe bilgi kutusu
        info_lines = []
        for res in vr:
            rn=res["run_no"]
            info_lines.append(
                f"Run {rn}:  n={res['n']}  R²={res['r2']:.4f}  "
                f"Slope={res['slope']:.4f}  Int={res['intercept']:.4f}"
            )
        info_txt = "\n".join(info_lines)
        ax.text(0.99, 0.03, info_txt,
            transform=ax.transAxes, fontsize=7.5,
            verticalalignment='bottom', horizontalalignment='right',
            color='#d0e0f0',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#1a3050',
                      edgecolor='#4a7090', alpha=0.9))

        ax.set_xlabel(self.T["x_logd"],color="#7090b0",fontsize=10)
        ax.set_ylabel(self.T["y_probit"],color="#7090b0",fontsize=10)
        ax.set_title(self.T["plot_title"],color="#FFC600",fontsize=12,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0",loc="upper left")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.pf)
        cv.draw(); cv.get_tk_widget().pack(fill="both",expand=True)

    def _dist(self):
        for w in self.df.winfo_children(): w.destroy()
        vr=[r for r in self.run_results if "error" not in r]
        if not vr: return
        stages=STAGE_ORDER+["MOC"]; x=np.arange(len(stages))
        n=len(vr); bw=0.8/max(n,1)
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        for i,res in enumerate(vr):
            c=CP[i%len(CP)]; ms=[res["masses"].get(s,0) for s in stages]
            off=(i-n/2+0.5)*bw
            ax.bar(x+off,ms,width=bw*0.88,color=c,alpha=0.85,label=f"Run {res['run_no']}")
        ax.set_xticks(x); ax.set_xticklabels(stages,rotation=30,ha="right",fontsize=9,color="#7090b0")
        ax.set_xlabel(self.T["x_stage"],color="#7090b0",fontsize=10)
        ax.set_ylabel(self.T["y_mass"],color="#7090b0",fontsize=10)
        ax.set_title(f"{self.T['dist_title']}  [{self.var_flow.get()} L/min]",
            color="#FFC600",fontsize=12,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
        ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.df)
        cv.draw(); cv.get_tk_widget().pack(fill="both",expand=True)

    def _sum(self):
        for w in self.sf.winfo_children(): w.destroy()
        vr=[r for r in self.run_results if "error" not in r]
        if not vr: return
        flow=int(self.var_flow.get())
        ctk.CTkLabel(self.sf,
            text=f"  {self.T['summary_hdr']}  |  {flow} L/min",
            font=ctk.CTkFont(size=13,weight="bold"),
            text_color="#FFC600",anchor="w").pack(fill="x",padx=6,pady=(8,4))
        params=[("metered",self.T["metered"],"mg"),
                ("delivered",self.T["delivered"],"mg"),
                ("fpd",self.T["fp_dose"],"mg"),
                ("fpf",self.T["fp_frac"],"%"),
                ("mmad",self.T["mmad"],"\u00b5m"),
                ("gsd",self.T["gsd"],""),
                ("slope",self.T["slope"],""),
                ("intercept",self.T["intercept"],""),
                ("r2",self.T["r2"],""),
                ("n",self.T["n_pts"],"")]
        hf=ctk.CTkFrame(self.sf,fg_color="#001a40",corner_radius=4)
        hf.pack(fill="x",padx=6,pady=1)
        ctk.CTkLabel(hf,text="Parameter",width=210,anchor="w",
            font=ctk.CTkFont(size=10,weight="bold"),
            text_color="#FFC600").pack(side="left",padx=8)
        for r in vr:
            ctk.CTkLabel(hf,text=f"Run {r['run_no']}",width=88,anchor="center",
                font=ctk.CTkFont(size=10,weight="bold"),
                text_color=CP[(r["run_no"]-1)%len(CP)]).pack(side="left",padx=2)
        for lbl in [self.T["mean"],self.T["sd"],self.T["rsd"]]:
            ctk.CTkLabel(hf,text=lbl,width=88,anchor="center",
                font=ctk.CTkFont(size=10,weight="bold"),
                text_color="#d0e0f0").pack(side="left",padx=2)
        sm=calc_summary(vr); hlk={"fpd","fpf","mmad","gsd"}
        for j,(key,lbl,unit) in enumerate(params):
            fr=ctk.CTkFrame(self.sf,
                fg_color="#090c14" if j%2==0 else "#141824",corner_radius=3)
            fr.pack(fill="x",padx=6,pady=1)
            ctk.CTkLabel(fr,text=f"{lbl} ({unit})" if unit else lbl,
                width=210,anchor="w",font=ctk.CTkFont(size=10),
                text_color="#7a9abf").pack(side="left",padx=8,pady=4)
            for r in vr:
                v=r.get(key,0); fmt=f"{float(v):.4f}" if isinstance(v,float) else str(v)
                ctk.CTkLabel(fr,text=fmt,width=88,anchor="center",
                    font=ctk.CTkFont(size=10),
                    text_color="#FFC600" if key in hlk else "#d0e0f0").pack(side="left",padx=2)
            m_,s_,rsd_=sm.get(key,(None,None,None))
            for sv in [m_,s_,rsd_]:
                t=f"{sv:.4f}" if sv is not None else "\u2014"
                ctk.CTkLabel(fr,text=t,width=88,anchor="center",
                    font=ctk.CTkFont(size=10),
                    text_color="#FFC600" if key in hlk else "#f0f0f0").pack(side="left",padx=2)

    def _pdf(self):
        if not self.run_results:
            messagebox.showwarning("",self.T["err_norun"]); return
        path=filedialog.asksaveasfilename(defaultextension=".pdf",
            filetypes=[("PDF","*.pdf")],
            initialfile=f"NGI_Report_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf")
        if not path: return
        flow=int(self.var_flow.get())
        meta={k:v.get() for k,v in self.mv.items()}
        try:
            make_pdf(path,self.run_results,meta,flow,self.T)
            self.lbl_status.configure(
                text=f"{self.T['pdf_saved']} {os.path.basename(path)}")
            if messagebox.askyesno("PDF","PDF a\u00e7\u0131ls\u0131n m\u0131?"):
                if os.name=="nt": os.startfile(path)
                else: os.system(f"xdg-open '{path}'")
        except Exception as e:
            messagebox.showerror("PDF Error", str(e))

    def _clr(self):
        for rw in self.rw:
            for v in rw["entries"].values(): v.set("0.000")
        self.run_results=[]
        for f in [self.rf,self.pf,self.df,self.sf]:
            for w in f.winfo_children(): w.destroy()
        self.lbl_status.configure(text=self.T["status_ready"])

    def _lang(self):
        self.lang="EN" if self.lang=="TR" else "TR"; self.T=L[self.lang]
        self.lbl_title.configure(text=self.T["title"])
        self.lbl_sub.configure(text=self.T["subtitle"])
        self.btn_lang.configure(text=self.T["lang_btn"])
        self.btn_add.configure(text=self.T["add_run"])
        self.btn_del.configure(text=self.T["del_run"])
        self.btn_calc.configure(text=self.T["calculate"])
        self.btn_clear.configure(text=self.T["clear"])
        self.btn_pdf.configure(text=self.T["export_pdf"])
        self.lbl_vr.configure(text=self.T["valid_range"])
        self.lbl_nb.configure(text=self.T["normalize_by"])
        self.rb1.configure(text=self.T["norm_ism"])
        self.rb2.configure(text=self.T["norm_total"])
        self._co()
        self.lbl_status.configure(text=self.T["status_ready"])

if __name__=="__main__":
    app=NGIApp(); app.mainloop()
