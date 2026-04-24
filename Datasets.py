import bone
pd = bone.pd
re = bone.re
hu = bone.hu
np = bone.np

asciiNorm = bone.asciiNorm

def plotTitleBar(cval, atypes, params):
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    w,h = (5, 0.8)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']
    if 'cval' in params:
        cval = params['cval']

    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig = bone.plt.figure(figsize=(w,h), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)
    nAt = len(cval[0])
    extent = [0, nAt, 0, 5]
    ax.axis(extent)
    cmap = bone.colors.ListedColormap(color_sch1)
    boundaries = range(len(color_sch1) + 1)
    norm = bone.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    #ax.imshow(cval, interpolation='none', cmap=cmap, \
    #                  norm=norm, extent=extent, aspect="auto")
    y = [0, 5]
    x = bone.np.arange(nAt + 1)
    ax.pcolormesh(x, y, cval, cmap=cmap, norm=norm, zorder=-1.0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(top=False, left=False, bottom=False, right=False)
    ax.set_xticks(bone.np.arange(0, nAt, 1))
    ax.grid(which='major', alpha=0.2, linestyle='-', linewidth=0.5,
            color='black', zorder=2.0)
    for edge, spine in ax.spines.items():
                spine.set_visible(False)
    divider = bone.make_axes_locatable(ax)
    width = bone.axes_size.AxesX(ax, aspect=1./20)
    spaceAnn = 70
    widthAnn = 3
    tAnn = 1
    if 'spaceAnn' in params:
        spaceAnn = params['spaceAnn']
    if 'widthAnn' in params:
        widthAnn = params['widthAnn']
    if 'tAnn' in params:
        tAnn = params['tAnn']
    pad = bone.axes_size.Fraction(0.1, width)
    lax = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
    lax.axison = False
    lax.axis(extent)
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.grid(False)
    lax.tick_params(top=False, left=False, bottom=False, right=False)
    if 'atypes' in params:
        atypes = params['atypes']
    bone.barTop(lax, atypes, color_sch1, params)
    return ax

bone.plotTitleBar = plotTitleBar

def plotTitleBarSingle(ana, ax=None, desc=None):
    expr = ana.f_ranks
    v1 = [expr[i-2] for i in ana.state[0]]
    v2 = [expr[i-2] for i in ana.state[1]]
    t, p = bone.ttest_ind(v1,v2, equal_var=False)
    expr = [expr[i-2] for i in ana.order]
    ana.cval = [[ana.aval[ana.order[i]] for i in bone.np.argsort(expr)]]
    ax = ana.printTitleBar({'w':4, 'h':0.5, 'spaceAnn': len(ana.order)/len(ana.atypes),
                            'tAnn': 1, 'widthAnn':1, 'ax':ax})
    if desc == None:
        desc = ana.h.getSource() + f"({len(ana.state[0])},{len(ana.state[1])})"
    ax.text(-1, 2, desc, horizontalalignment='right', verticalalignment='center')
    actual = ana.cval[0]
    score = [expr[i] for i in bone.np.argsort(expr)]
    fpr, tpr, thrs = bone.roc_curve(actual, score, pos_label=1)
    auc = bone.auc(fpr, tpr)
    roc_auc = "%.2f %.3g" % (auc, p)
    print(roc_auc)
    ax.text(len(ana.cval[0]), 2, roc_auc)
    return [auc, p, t]
bone.plotTitleBarSingle = plotTitleBarSingle

def adjustFigure(ana, fig, xlabel="Score", wfactor=1):
    rocauc = ana.getROCAUC()
    ax = fig.axes[1]
    #ax.set_xlim([-5, 10])
    ax.set_title(f"ROC-AUC = {rocauc}", fontsize=16)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=20)
    #fig.patch.set_facecolor('#FFF7CD')
    #ax.set_facecolor('#FFF7CD')
    for tobj in ax.texts:
        tobj.set_fontsize(16)
        #tobj.set_color('black')
    #for tobj in fig.findobj(plt.Text):
    #    tobj.set_color('black')        
    fw, fh = fig.get_size_inches()
    bbox = fig.axes[2].get_position()
    aw, ah = bbox.width * fw, bbox.height * fh
    xlim = fig.axes[2].get_xlim()
    ylim = fig.axes[2].get_ylim()
    for p in fig.axes[2].patches:
        bbox1 = p.get_bbox()
        pw = bbox1.width/(xlim[1]-xlim[0]) * aw * 72
        ph = bbox1.height/(ylim[1]-ylim[0]) * ah * 72
        sw = ph/72/aw*(xlim[1]-xlim[0])/wfactor
        p.set_width(sw)
    return
bone.adjustFigure = adjustFigure

def getMSigDB(gs):
    url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=" + gs + "&fileType=txt"
    df = pd.read_csv(url, sep="\t")
    df.columns.values[0] = 'ID'
    l1 = [list(df.ID[1:])]
    wt1 = [1]
    return wt1, l1
bone.getMSigDB = getMSigDB

def prepareData(self, dbid, cfile = "data/explore.conf"):
    self.db = hu.Database(cfile)
    self.dbid = dbid
    self.h = hu.Hegemon(self.db.getDataset(self.dbid))
    self.h.init()
    self.h.initPlatform()
    self.h.initSurv()
    self.num = self.h.getNum()
    self.start = self.h.getStart()
    self.end = self.h.getEnd()
    self.name = self.h.rdataset.getName()
    self.source = self.h.getSource()
    self.headers = self.h.headers
    self.hhash = {}
    for i in range(len(self.headers)):
        self.hhash[self.headers[i]] = i
    return
bone.IBDAnalysis.prepareData = prepareData

def getSMheart(self, tn=1):
    self.prepareData("SM1", cfile="data/explore.conf")
    atype = self.h.getSurvName("c Title")
    atypes = ['WT', 'KO', 'KO+CST']
    ahash = {'W2_S2':0, 'W3_S3':0, 'W1_S1':0,
            'C1_S10':1, 'C2_S11':1, 'C3_S12':1, 'C5_S14':1,
            'CC1_S15':2, 'CC2_S16':2, 'CC4_S18':2, 'CC5_S19':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSMheart = getSMheart

def getMorley2019(self, tn=1):
    self.prepareData("COV109")
    atype = self.h.getSurvName("c etiology")
    atypes = ['N', 'DCM', 'HCM', 'PPCM']
    ahash = {'Dilated cardiomyopathy (DCM)':1,
            'Non-Failing Donor':0,
            'Hypertrophic cardiomyopathy (HCM)':2,
            'Peripartum cardiomyopathy (PPCM)':3}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['H', 'CM']
        ahash = {'Dilated cardiomyopathy (DCM)':1,
                'Non-Failing Donor':0,
                'Hypertrophic cardiomyopathy (HCM)':1,
                'Peripartum cardiomyopathy (PPCM)':1}
    if (tn == 3):
        atype = self.h.getSurvName("c race")
        atypes = ['W', 'B']
        ahash = {'Caucasian':0, 'African American':1}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorley2019 = getMorley2019

def getHannenhalli2006(self, tn=1):
    self.prepareData("COV113")
    atype = self.h.getSurvName("c src1")
    atypes = ['NF', 'LV']
    ahash = {'explanted heart tissue at time of cardiac transplantation':1,
             'unused donor heart with normal LV function':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHannenhalli2006 = getHannenhalli2006

def getDorsey2016(self, tn=1):
    self.prepareData("MUSCLE23")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'DMD']
    ahash = {'control_muscle biopsy':0, 'DMD_muscle biopsy':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDorsey2016 = getDorsey2016

def getDeMasi2024Mm(self, tn=1):
    self.prepareData("MUSCLE24")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_rep.", "", str(k)) for k in atype]
    atypes = ['TAC', 'TACH', 'HC', 'HCH', 'GC', 'GCH', 'O']
    ahash = {'gastrocnemius_BL10_ctrl':4, 'gastrocnemius_mdx_ctrl':4,
            'gastrocnemius_mdx_CHP':5, 'heart_BL10_ctrl':2,
            'heart_mdx_ctrl':2, 'heart_mdx_CHP':3, '':6,
            'tibialis_aged_CHP':1, 'tibialis_aged_ctrl':0}
    if tn == 2:
        atypes = ['BL10', 'DMD']
        ahash = {'heart_BL10_ctrl':0, 'heart_mdx_ctrl':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeMasi2024Mm = getDeMasi2024Mm

def getMarini2022(self, tn=1):
    self.prepareData("CA10")
    atype = self.h.getSurvName("c genotype/variation")
    atypes = ['H', 'DMD']
    ahash = {'Contains specific DMD point mutation':1,
            'DMD point mutation was corrected via CRISPR/Cas technology':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMarini2022 = getMarini2022

def getGalindo2016(self, tn=1):
    self.prepareData("CA11")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" biolo.*", "", str(k)) for k in atype]
    atypes = ['VNY', 'VNM', 'VGY', 'VGM', 'VGO',
            'SNY', 'SNM', 'SGY', 'SGM', 'SGO']
    ahash = {'left ventricle   Normal   young':0,
            'left ventricle   Normal   young ':0,
            'left ventricle   Normal   middle':1,
            'left ventricle   GRMD   young':2,
            'left ventricle   GRMD   middle':3,
            'left ventricle   GRMD   old':4,
            'skeletal muscle   Normal   young':5,
            'skeletal muscle   Normal   middle':6,
            'skeletal muscle   GRMD   young':7,
            'skeletal muscle   GRMD   middle':8,
            'skeletal muscle   GRMD   old':9}
    if tn == 2:
        atypes = ['SNY', 'SGY']
        ahash = {'skeletal muscle   Normal   young':0,
                'skeletal muscle   GRMD   young':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGalindo2016 = getGalindo2016

def getLi2020(self, tn=1):
    self.prepareData("COV107")
    atype = self.h.getSurvName("c src1")
    atypes = ['LV', 'ARV']
    ahash = {"ARVC patients' heart LV":0, 'heart RV tissue':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2020 = getLi2020

def getGlobal(self, tn=1):
    self.prepareData("GL1")
    atype = self.h.getSurvName('c TissueState')
    atypes = ['nan', 'solid', 'cellline', 'liquid', 'bonemarrow',
            'mix', 'unknown', 'lymphoma', 'lymph',
            'cellline_p3', 'saliva']
    ahash = {'':0}
    if tn == 2:
        atypes = ['solid', 'liquid', 'cellline']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobal = getGlobal

def getKatoh2025Heart(self, tn=1, ta=0):
    self.prepareData("CA48")
    atype = self.h.getSurvName('c Type')
    atypes = ['CTRL', 'DCM']
    ahash = {}
    if (tn == 2): # Select Only males
        gsmid = self.h.getSurvName('c gsmid')
        atype = [None if gsmid[k] == 'GSM9102209' else atype[k]
                for k in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKatoh2025Heart = getKatoh2025Heart

def getRen2020Mm(self, tn=1):
    self.prepareData("COV108.2")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'C', 'CM']
    ahash = {'TAC sham CM':1, 'TAC2W CM':2, 'TAC5W CM':2,
            'TAC8W CM':2, 'TAC11W CM':2, 'normal heart CM':0}
    if (tn == 2):
        atypes = ['C', 'CM']
        ahash = {'TAC sham CM':0, 'TAC2W CM':1, 'TAC5W CM':1,
                'TAC8W CM':1, 'TAC11W CM':1, 'normal heart CM':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRen2020Mm = getRen2020Mm

def getGrueter2025heartMm(self, tn=1, ta=0):
    self.prepareData("CA50")
    atype = self.h.getSurvName('c Type')
    atypes = ['C', 'CM']
    ahash = {'Ventricle, Neg Cre':0, 'Ventricle, Pos Cre':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGrueter2025heartMm = getGrueter2025heartMm

def getPavelec2025heartMm(self, tn=1, ta=0):
    self.prepareData("CA51")
    atype = self.h.getSurvName('c Type')
    atypes = ['SW', 'SK', 'IW', 'IK', 'W', 'K']
    ahash = {'Cardiomyocyte-Saline-WT':0, 'Cardiomyocyte-Saline-KO':1,
            'Cardiomyocyte-ISO-WT':2, 'Cardiomyocyte-ISO-KO':3,
            'WholeHeart-ISO-KO':5, 'WholeHeart-ISO-WT':4}
    if (tn == 2):
        atypes = ['C', 'CM']
        ahash = {'Cardiomyocyte-Saline-WT':1, 'Cardiomyocyte-Saline-KO':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPavelec2025heartMm = getPavelec2025heartMm

def getBakshi2025heartMm(self, tn=1, ta=0):
    self.prepareData("CA54")
    atype = self.h.getSurvName('c Type')
    atypes = ['CSV', 'CTV', 'CTM', 'SSV', 'STV', 'STM']
    ahash = {'Whole heart tissue,Control-Sham-Vehicle':0,
            'Whole heart tissue,Control-TAC-Vehicle':1,
            'Whole heart tissue,Control-TAC-MitoQ':2,
            'Whole heart tissue,cSirt4-Tg-Sham-Vehicle':3,
            'Whole heart tissue,cSirt4-Tg-TAC-Vehicle':4,
            'Whole heart tissue,cSirt4-Tg-TAC-MitoQ':5}
    if (tn == 2):
        atypes = ['C', 'CM']
        ahash = {'Whole heart tissue,Control-Sham-Vehicle':0,
                'Whole heart tissue,Control-TAC-Vehicle':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBakshi2025heartMm = getBakshi2025heartMm

def getPeng2020htnRt(self, tn=1):
    self.prepareData("HRT6")
    atype = self.h.getSurvName("c group")
    atypes = ['W', 'S', 'S+Q']
    ahash = {'Hypertension':1, 'QDG Treatment':2, 'Normal':0}
    if tn == 2:
        atypes = ['W', 'S']
        ahash = {'Hypertension':1, 'Normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeng2020htnRt = getPeng2020htnRt

def getKhassafi2023htnRt(self, tn=1, ta=0):
    self.prepareData("CA40.2")
    atype = self.h.getSurvName('c condition')
    atypes = ['N', 'C', 'D', 'E', 'L']
    ahash = {'Control RV':0, 'control RV':0,
            'Compensated RV':1, 'Decompensated RV':2,
            'early decompensated':3, 'late decompensated':4}
    if tn == 2:
        atype = self.h.getSurvName('c tissue type')
        atypes = ['N', 'C', 'D']
        ahash = {'Control RV':0, 'Compensated RV':1, 'Decompensated RV':2}
    if tn == 3:
        atypes = ['N', 'C']
        ahash = {'Control RV':0, 'control RV':0, 'Compensated RV':1}
    if tn == 4:
        atypes = ['N', 'D']
        ahash = {'Control RV':0, 'control RV':0, 'Decompensated RV':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKhassafi2023htnRt = getKhassafi2023htnRt

