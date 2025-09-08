def compute_jet_id_single(jet):
    """Return (passTight, passTightLepVeto) for one jet (AK4 or FatJet).

    """
    eta = abs(getattr(jet, "eta", 0.0))
    chMult = int(getattr(jet, "chMultiplicity", 0))
    neMult = int(getattr(jet, "neMultiplicity", 0))
    neHEF  = float(getattr(jet, "neHEF",  0.0))
    neEmEF = float(getattr(jet, "neEmEF", 0.0))
    chHEF  = float(getattr(jet, "chHEF",  0.0))
    muEF   = float(getattr(jet, "muEF",   0.0))
    chEmEF = float(getattr(jet, "chEmEF", 0.0))

    passTight = (
        ((eta <= 2.6) and (neHEF < 0.99) and (neEmEF < 0.9) and ((chMult + neMult) > 1) and (chHEF > 0.01) and (chMult > 0))
        or
        ((2.6 < eta <= 2.7) and (neHEF < 0.90) and (neEmEF < 0.99))
        or
        ((2.7 < eta <= 3.0) and (neHEF < 0.99))
        or
        ((eta > 3.0) and (neMult >= 2) and (neEmEF < 0.4))
    )

    # LepVeto applies only for |eta|<=2.7; beyond that it's identical to passTight
    passTightLepVeto = (passTight and (muEF < 0.8) and (chEmEF < 0.8)) if (eta <= 2.7) else passTight
    return passTight, passTightLepVeto
