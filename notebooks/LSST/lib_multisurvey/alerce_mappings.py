# ALeRCE database queries
# Authors: Alejandra Muñoz Arancibia et al.
#
# Definitions for survey-related properties


sid_num2str = {
    0: "ZTF",
    1: "LSST diaObject",
    2: "LSST ssObject",
}

sid_survey = {
    0: "ZTF",
    1: "LSST",
    2: "LSST",
}

sid_map_cols = {
    0: {
        "oid": "oid",
        "candid": "candid",
        "ssid": "ssnamenr",
        "sscandid": "candid",
        "mjd": "mjd",
        "band": "band",
        "flux_diff": "fluxdiff_nJy",
        "flux_diff_err": "fluxdiff_err_nJy",
        "flux_sci": "fluxsci_nJy",
        "flux_sci_err": "fluxsci_err_nJy",
    },
    1: {
        #'oid': 'diaObjectId',
        #'candid': 'diaSourceId',
        #'mjd': 'midpointMjdTai',
        "oid": "oid",
        "candid": "measurement_id",
        "ssid": "ssObjectId",
        "sscandid": "ssSourceId",
        "mjd": "mjd",
        "band": "band_name",
        "flux_diff": "psfFlux",
        "flux_diff_err": "psfFluxErr",
        "flux_sci": "scienceFlux",
        "flux_sci_err": "scienceFluxErr",
    },
    2: {
        #'oid': 'diaObjectId',
        #'candid': 'diaSourceId',
        #'mjd': 'midpointMjdTai',
        "oid": "oid",
        "candid": "measurement_id",
        "ssid": "ssObjectId",
        "sscandid": "ssSourceId",
        "mjd": "mjd",
        "band": "band_name",
        "flux_diff": "psfFlux",
        "flux_diff_err": "psfFluxErr",
        "flux_sci": "scienceFlux",
        "flux_sci_err": "scienceFluxErr",
    },
}

sid_bands = {
    0: ["g", "r", "i"],
    1: ["u", "g", "r", "i", "z", "y"],
    2: ["u", "g", "r", "i", "z", "y"],
    
}

sid_map_bands = {
    0: {1: "g", 2: "r", 3: "i"},
    1: {6: "u", 1: "g", 2: "r", 3: "i", 4: "z", 5: "y"},
    2: {6: "u", 1: "g", 2: "r", 3: "i", 4: "z", 5: "y"},
    
}

sid_map_cols_forced2det = {
    0: {
        "mag": "magpsf",
        "e_mag": "sigmapsf",
        "mag_corr": "magpsf_corr",
        "e_mag_corr_ext": "sigmapsf_corr_ext",
    },
}

image_subtypes = [
    "Science",
    "Template",
    "Difference",
]

# Only for LSST
image_planes = ["flux", "variance", "mask"]

tables_epochs = [
    "detections",
    "forced_photometry",
    # "non_detections",
]
