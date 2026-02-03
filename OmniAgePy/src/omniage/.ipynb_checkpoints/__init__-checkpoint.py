__version__ = "0.1.0"

from .predict import cal_epimarker, list_available_clocks, get_clock_coefs,get_clock_coefs

from .models import (
    # --- Epigenetic (DNAm) ---
    Horvath2013, Hannum, Lin, VidalBralo, ZhangClock, Horvath2018, 
    Bernabeu_cAge, CorticalClock, PedBE, CentenarianClock_40, CentenarianClock_100,
    Retro_age_V1, Retro_age_V2,
    Zhang10, PhenoAge, DunedinPACE, GrimAge1, GrimAge2, DNAmFitAge, IC_Clock,
    EpiTOC1, EpiTOC2, EpiTOC3, StemTOC, StemTOCvitro, RepliTali, HypoClock, EpiCMIT_Hypo,
    EpiCMIT_Hyper,DNAmTL, PCHorvath2013, PCHorvath2018, PCHannum, PCPhenoAge, PCDNAmTL, PCGrimAge1,
    ABEC, eABEC, cABEC, PipekElasticNet, PipekFilteredh, PipekRetrainedh, WuClock,
    # --- Causal / Stochastic ---
    CausalAge, DamAge, AdaptAge, StocH, StocP, StocZ,
    
    # --- Systems Age ---
    SystemsAge,
    
    # --- Surrogate Biomarkers ---
    CompCRP, CompCHIP, EpiScores, CompIL6,
    McCartneyBMI, McCartneySmoking, McCartneyAlcohol, McCartneyEducation,
    McCartneyTotalCholesterol, McCartneyHDL, McCartneyLDL,
    McCartneyTotalHDLRatio, McCartneyWHR, McCartneyBodyFat,
    
    # --- Transcriptomic ---
    scImmuAging, BrainCTClock, PASTA_Clock,
    
    # --- CTS & Pan-Mammalian & Ensemble ---
    NeuIn, GliaIn, Brain, NeuSin, GliaSin, Hep, Liver,
    PanMammalianUniversal, PanMammalianBlood, PanMammalianSkin,
    EnsembleAgeHumanMouse,EnsembleAgeStatic,EnsembleAgeDynamic,
    
    # --- Gestational & CTF ---
    BohlinGA, EPICGA, LeeControl,LeeRobust,LeeRefinedRobust, KnightGA, MayneGA,
    DNAmCTFClock,

    # --- Disease Risk ---
    CompSmokeIndex
)



# 3. define __all__
__all__ = [
    "cal_epimarker",
    
    # DNAm Clocks
    "Horvath2013", "Hannum", "Lin", "VidalBralo", "ZhangClock", "Horvath2018", 
    "Bernabeu_cAge", "CorticalClock", "PedBE", "CentenarianClock_40", "CentenarianClock_100",
    "Retro_age_V1", "Retro_age_V2",
    "ABEC","eABEC","cABEC","PipekElasticNet","PipekFilteredh","PipekRetrainedh","WuClock",
    "Zhang10", "PhenoAge", "DunedinPACE", "GrimAge1", "GrimAge2", "DNAmFitAge", "IC_Clock",
    "EpiTOC1", "EpiTOC2", "EpiTOC3", "StemTOC", "StemTOCvitro", "RepliTali", "HypoClock", "EpiCMIT_Hypo",
    "EpiCMIT_Hyper","DNAmTL","PCHorvath2013", "PCHorvath2018", "PCHannum", "PCPhenoAge", "PCDNAmTL", "PCGrimAge1",
    
    # Advanced
    "CausalAge", "DamAge", "AdaptAge", "StocH", "StocP", "StocZ", "SystemsAge",
    
    # Surrogates
    "CompCRP", "CompCHIP", "EpiScores", "CompIL6",
    "McCartneyBMI", "McCartneySmoking", "McCartneyAlcohol", "McCartneyEducation",
    "McCartneyTotalCholesterol", "McCartneyHDL", "McCartneyLDL",
    "McCartneyTotalHDLRatio", "McCartneyWHR", "McCartneyBodyFat",
    
    # Transcriptomic
    "scImmuAging", "BrainCTClock", "PASTA_Clock",
    
    # Others
    "NeuIn", "GliaIn", "Brain", "NeuSin", "GliaSin", "Hep", "Liver",
    "PanMammalianUniversal", "PanMammalianBlood", "PanMammalianSkin",
    "EnsembleAgeHumanMouse","EnsembleAgeStatic","EnsembleAgeDynamic",
    "BohlinGA", "EPICGA", "LeeControl","LeeRobust","LeeRefinedRobust", "KnightGA", "MayneGA",
    "DNAmCTFClock","CompSmokeIndex"
]