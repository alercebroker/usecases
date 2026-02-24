# ALeRCE multisurvey plotting toolbox
# Authors: Alejandra Muñoz Arancibia et al.
#
# Based on many notebooks by the ALeRCE collaboration
#
# Functions for plotting light curves and stamps


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.wcs import WCS
import sys

import alerce_mappings as am


client_path = "/home/ammunoz1/Desktop/Astro2025/ALeRCE/alercebroker/alerce_client/"
sys.path.append(client_path)
from alerce.core import Alerce

alerce_client = Alerce()


lc_colors = {
    0: {
        "g": "#56E03A",
        "r": "#D42F4B",
        "i": "gold",
    },
    ## Based on https://github.com/lsst/utils/commit/8821c3fa7c1f637ec19bab1ca06eacbf6db81cfc#diff-88ddc631311241ce060988ae76aa32db53575fec9597e162f70bb6f8fd975c57
    # 1: {
    #    "u": "#0c71ff",
    #    "g": "#49be61",
    #    "r": "#c61c00",
    #    "i": "#ffc200",
    #    "z": "#f341a2",
    #    "y": "#5d0000",
    # },
    # Used in ALeRCE Explorer
    1: {
        "u": "#56B4E9",  # sky blue
        "g": "#009E73",  # bluish green
        "r": "#D55E00",  # vermillion
        "i": "#E69F00",  # orange
        "z": "#CC79A7",  # reddish purple
        "y": "#0072B2",  # blue
    },
    2: {
        "u": "#56B4E9",  # sky blue
        "g": "#009E73",  # bluish green
        "r": "#D55E00",  # vermillion
        "i": "#E69F00",  # orange
        "z": "#CC79A7",  # reddish purple
        "y": "#0072B2",  # blue
    },
}

lc_markers = {
    0: {
        "g": "o",
        "r": "s",
        "i": "X",
    },
    1: {
        "u": "X",
        "g": "D",
        "r": "p",
        "i": "o",
        "z": "s",
        "y": "*",
    },
    2: {
        "u": "X",
        "g": "D",
        "r": "p",
        "i": "o",
        "z": "s",
        "y": "*",
    },
}
## Based on https://github.com/lsst/utils/blob/main/python/lsst/utils/plotting/figures.py
# lc_markers = {
#    1: {
#        'u': 'o',
#        'g': '^',
#        'r': 'v',
#        'i': 's',
#        'z': '*',
#        'y': 'p',
#    }
# }

lc_sizes = {
    0: {
        "g": 40,
        "r": 40,
        "i": 40,
    },
    1: {
        "u": 50,
        "g": 40,
        "r": 60,
        "i": 50,
        "z": 40,
        "y": 70,
    },
    2: {
        "u": 50,
        "g": 40,
        "r": 60,
        "i": 50,
        "z": 40,
        "y": 70,
    },
}
# lc_sizes = {
#    1: {
#        'u': 60,
#        'g': 60,
#        'r': 60,
#        'i': 60,
#        'z': 60,
#        'y': 60,
#    }
# }
lc_sizes_forced = {
    0: {key: val - 20 for key, val in lc_sizes[0].items()},
    1: {key: val - 20 for key, val in lc_sizes[1].items()},
    2: {key: val - 20 for key, val in lc_sizes[2].items()},
}

# Jan 1 from 2023 to 2027
mjd_jan1 = [59945.0, 60310.0, 60676.0, 61041.0, 61406.0]
year_i = 2023
year_f = 2027
yr = pd.DataFrame(mjd_jan1, columns=["mjd_jan1"])
yr["label"] = np.array(range(year_i, year_f + 1)).astype(str)


def plot_stamps_fromapi(df_objs=None, survey="lsst", sid=1):
    df = pd.DataFrame()

    col_id = am.sid_map_cols[sid]["oid"]
    col_candid = am.sid_map_cols[sid]["candid"]

    for i, (oid, candid) in enumerate(zip(df_objs.index.tolist(), df_objs[col_candid])):
        # print(i, oid, candid)

        params = {
            "survey": survey,
            col_id: str(oid),
            col_candid: str(candid),
        }
        # print(params)

        stamps = alerce_client.plot_stamps(**params)


def plot_stamps(
    row,
    imtypes=am.image_subtypes,
    plane="flux",
    title="",
    use_wcs=False,
    figsize=[5.5, 5.5],
    fontsize=None,
    title_offset=0.72,
    namefig=None,
):
    if fontsize is not None:
        plt.rcParams.update({"font.size": fontsize})

    ncols = len(imtypes)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrows=1, ncols=ncols, figure=fig)

    fig.suptitle(title, y=title_offset)

    for i, imtype in enumerate(imtypes):
        im = row[plane + "_" + imtype + "_data"]
        if plane == "flux":
            im = np.arcsinh(im)

        if use_wcs and plane == "flux":
            w = WCS(row[plane + "_" + imtype + "_header"])
            ax = fig.add_subplot(gs[int(i / ncols), i % ncols], projection=w)
        else:
            ax = fig.add_subplot(gs[int(i / ncols), i % ncols])

        ax.imshow(im, cmap="viridis")
        ax.set_title(imtype)

        if use_wcs and plane == "flux":
            lon = ax.coords[0]
            lat = ax.coords[1]
            lon.set_axislabel("RA [deg]")
            lat.set_axislabel("Dec [deg]")
            # lon.set_ticks_visible(False)
            # lon.set_ticklabel_visible(False)
            # lat.set_ticks_visible(False)
            # lat.set_ticklabel_visible(False)
        else:
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)

    plt.tight_layout()

    if namefig is not None:
        fig.savefig(namefig, bbox_inches="tight")

    if fontsize is not None:
        plt.rcParams.update({"font.size": 12})


def query_stamps(
    df_objs=None,
    survey="lsst",
    sid=1,
    imtypes=am.image_subtypes,
    include_variance_and_mask=False
):
    df = pd.DataFrame()

    col_id = am.sid_map_cols[sid]["oid"]
    col_candid = am.sid_map_cols[sid]["candid"]

    for i, (oid, candid) in enumerate(zip(df_objs.index.tolist(),
                                          df_objs[col_candid])):
        # print(i, oid, candid)

        params = {
            "survey": survey,
            col_id: oid,
            "candid": candid,
        }
        # print(params)

        try:
            stamps = alerce_client.get_stamps(
                **params, include_variance_and_mask=include_variance_and_mask
            )
            # print(stamps)
        except:
            stamps = None

        if stamps is None:
            continue

        obj = {col_id: oid, col_candid: candid}
        # print(stamps)

        planes = (
            am.image_planes
            if isinstance(stamps["cutout" + imtypes[0]], list)
            else ["flux"]
        )

        for j, plane in enumerate(planes):
            for imtype in imtypes:
                data = (
                    stamps["cutout" + imtype][j]
                    if isinstance(stamps["cutout" + imtypes[0]], list)
                    else stamps["cutout" + imtype]
                )

                obj[plane + "_" + imtype + "_data"] = [data.data]
                obj[plane + "_" + imtype + "_header"] = [data.header]

        obj = pd.DataFrame(obj)
        df = pd.concat([df, obj], axis=0)
        del obj

    df.set_index(col_id, inplace=True)

    return df


# Based on https://github.com/alercebroker/pipeline/blob/main/lc_classifier/lc_classifier/utils.py
def flux_err_2_mag_err(flux_err, flux):
    return (2.5 * flux_err) / (np.log(10.0) * flux)


# Based on https://github.com/lsst/tutorial-notebooks/blob/main/DP0.2/01_Introduction_to_DP02.ipynb
def fluxnjy2mag(flux):
    return 31.4 - 2.5 * np.log10(flux)


def absmag(appmag, distmod=None):
    return appmag - distmod


def mag2fluxnjy(mag):
    return 10 ** (-(mag - 31.4) / 2.5)


def mag_err_2_fluxnjy_err(mag_err, mag):
    return np.log(10.0) * mag2fluxnjy(mag) / 2.5 * mag_err


def limit_epochs(d_objs=None, mjd_lims=None):
    d_objs = d_objs.copy()

    for sid in d_objs.keys():
        col_mjd_this = am.sid_map_cols[sid]["mjd"]

        for oid in d_objs[sid].keys():
            for level in ["detections", "forced_photometry", "non_detections"]:
                if len(d_objs[sid][oid][level]) > 0:
                    mask = (d_objs[sid][oid][level][col_mjd_this] >= mjd_lims[0]) \
                        & (d_objs[sid][oid][level][col_mjd_this] <= mjd_lims[1])
                    d_objs[sid][oid][level] = d_objs[sid][oid][level][mask].copy()

    return d_objs


def group_data_sid_oid(
    d_objs={},
    sid=None,
    oid=None,
    df_dets=pd.DataFrame(),
    df_forced=pd.DataFrame(),
    df_nondets=pd.DataFrame(),
    df_features=pd.DataFrame(),
    dict_others={},
    lc_kwargs={},
):
    d_objs = d_objs.copy()

    if sid not in d_objs.keys():
        d_objs[sid] = {}

    if oid not in d_objs[sid].keys():
        d_objs[sid][oid] = {}

    d_objs[sid][oid]["detections"] = df_dets.copy()
    d_objs[sid][oid]["forced_photometry"] = df_forced.copy()
    d_objs[sid][oid]["non_detections"] = df_nondets.copy()
    d_objs[sid][oid][
        "features"
    ] = df_features.copy()  # Not used yet - TODO choose between pivot or not
    d_objs[sid][oid]["others"] = dict_others.copy()  # Not used yet
    d_objs[sid][oid]["lc_kwargs"] = lc_kwargs.copy()

    return d_objs


def init_lc_kwargs():
    lc_kwargs = {}

    keys = [
        "show_forced",
        "show_flux",
        "show_diff",
        "show_sci",
    ]
    for key in keys:
        lc_kwargs[key] = True

    keys = [
        "show_absolute_mag",
        "show_folded",
    ]
    for key in keys:
        lc_kwargs[key] = False

    keys = ["z_obj", "folded_period"]
    for key in keys:
        lc_kwargs[key] = np.nan
    
    keys = ["cosmo"]
    for key in keys:
        lc_kwargs[key] = None

    lc_kwargs["title_exts_flux_diff"] = ", difference flux"
    lc_kwargs["title_exts_flux_sci"] = ", science flux"
    lc_kwargs["title_exts_mag_diff"] = ", difference magnitude"
    lc_kwargs["title_exts_mag_sci"] = ", apparent magnitude"
    lc_kwargs["title_exts_absmag_diff"] = ", absolute magnitude (from diff)"
    lc_kwargs["title_exts_absmag_sci"] = ", absolute magnitude"

    return lc_kwargs


# TODO: add non_detections in plot
# TODO: add low-S/N forced photometry epochs as upper limits in plot
# TODO: check if absolute magnitude plots work
def plot_lc(
    d_objs={}, mjd_lims=None, y_lims=None, lc_params={}, title_exts="",
    namefig=None
):
    lc_params_default = {
        "from_tap": False,
        "use_mag": False,
        "use_flux": True,
        "show_yr": False,
        "figwidth": 9,
        "pheight": 3.5,
        "fontsize": None,
        "alpha": 0.7,
        "alpha_forced": 0.5,
        "offset_mjd": 5.0,
    }
    lc_params = lc_params.copy()
    for param in lc_params_default.keys():
        if param not in lc_params.keys():
            lc_params[param] = lc_params_default[param]

    from_tap = lc_params["from_tap"]
    use_mag = lc_params["use_mag"]
    use_flux = lc_params["use_flux"]
    show_yr = lc_params["show_yr"]
    figwidth = lc_params["figwidth"]
    pheight = lc_params["pheight"]
    fontsize = lc_params["fontsize"]
    alpha = lc_params["alpha"]
    alpha_forced = lc_params["alpha_forced"]
    offset_mjd = lc_params["offset_mjd"]

    if fontsize is not None:
        plt.rcParams.update({"font.size": fontsize})

    d_objs = d_objs.copy()

    if mjd_lims is not None:
        d_objs = limit_epochs(d_objs=d_objs, mjd_lims=mjd_lims)

    nrows = 0

    level = "detections"
    for sid in d_objs.keys():
        for oid in d_objs[sid].keys():
            if len(d_objs[sid][oid][level]) > 0:
                if d_objs[sid][oid]["lc_kwargs"]["show_diff"]:
                    nrows += 1
                if d_objs[sid][oid]["lc_kwargs"]["show_sci"]:
                    nrows += 1

    fig = plt.figure(figsize=(figwidth, pheight * nrows))
    gs = GridSpec(nrows=nrows, ncols=1, figure=fig)

    i = 0

    for sid in d_objs.keys():
        bands = am.sid_bands[sid]
        col_band = am.sid_map_cols[sid]["band"]
        col_mjd = am.sid_map_cols[sid]["mjd"]
        colors = lc_colors[sid]
        markers = lc_markers[sid]
        sizes = lc_sizes[sid]
        sizes_forced = lc_sizes_forced[sid]

        for oid in d_objs[sid].keys():
            title = am.sid_num2str[sid] + ", " + am.sid_map_cols[sid]["oid"] \
                    + "=" + str(oid)

            df_dets = d_objs[sid][oid]["detections"].copy()

            if len(df_dets) == 0:
                continue

            use_flux = d_objs[sid][oid]["lc_kwargs"]["show_flux"]
            use_absmag = d_objs[sid][oid]["lc_kwargs"]["show_absolute_mag"]
            use_folded = d_objs[sid][oid]["lc_kwargs"]["show_folded"]

            z_obj = d_objs[sid][oid]["lc_kwargs"]["z_obj"]
            cosmo = d_objs[sid][oid]["lc_kwargs"]["cosmo"]
            period = d_objs[sid][oid]["lc_kwargs"]["folded_period"]

            if use_absmag and ~np.isnan(z_obj):
                distmod = cosmo.distmod(z_obj).value

            for light_type in ["diff", "sci"]:
                show_light_type = d_objs[sid][oid]["lc_kwargs"]["show_" + light_type]

                if not show_light_type:
                    continue

                ax_i = plt.subplot(gs[i])

                col_y = am.sid_map_cols[sid]["flux_" + light_type]
                col_yerr = am.sid_map_cols[sid]["flux_" + light_type + "_err"]
                
                if sid != 0 and from_tap:
                    col_y = col_y.lower()
                    col_yerr = col_yerr.lower()

                for band in bands:
                    mask = df_dets[col_band] == band
                    if len(df_dets[mask]) > 0:
                        if use_flux:
                            y = df_dets[col_y][mask]
                            yerr = df_dets[col_yerr][mask]
                        else:
                            y = fluxnjy2mag(df_dets[col_y][mask])
                            yerr = flux_err_2_mag_err(
                                df_dets[col_yerr][mask], df_dets[col_y][mask].abs()
                            )

                            if use_absmag and z_obj is not None:
                                y = absmag(y, distmod=distmod)

                        if use_folded and period is not None:
                            phase = np.mod(df_dets[col_mjd][mask], period) / period

                            for nphase in [0, 1]:
                                if nphase == 0:
                                    label = "%s" % band
                                else:
                                    label = None

                                ax_i.errorbar(
                                    phase + nphase,
                                    y,
                                    yerr=yerr,
                                    alpha=alpha,
                                    c=colors[band],
                                    linestyle="None",
                                )
                                ax_i.scatter(
                                    phase + nphase,
                                    y,
                                    s=sizes[band],
                                    alpha=alpha,
                                    c=colors[band],
                                    marker=markers[band],
                                    label=label,
                                    linestyle="None",
                                )
                        else:
                            ax_i.errorbar(
                                df_dets[col_mjd][mask],
                                y,
                                yerr=yerr,
                                alpha=alpha,
                                c=colors[band],
                                linestyle="None",
                            )
                            ax_i.scatter(
                                df_dets[col_mjd][mask],
                                y,
                                s=sizes[band],
                                alpha=alpha,
                                c=colors[band],
                                marker=markers[band],
                                label="%s" % band,
                                linestyle="None",
                            )

                if mjd_lims is None and len(df_dets[df_dets[col_y].notna()]) == 1:
                    ax_i.set_xlim(
                        [
                            df_dets.iloc[0][col_mjd] - offset_mjd,
                            df_dets.iloc[0][col_mjd] + offset_mjd,
                        ]
                    )

                df_forced = d_objs[sid][oid]["forced_photometry"].copy()
                use_forced = d_objs[sid][oid]["lc_kwargs"]["show_forced"]

                if len(df_forced) > 0 and use_forced:
                    # display(df_forced)
                    for band in bands:
                        mask = df_forced[col_band] == band
                        if len(df_forced[mask]) > 0:
                            if use_flux:
                                y = df_forced[col_y][mask]
                                yerr = df_forced[col_yerr][mask]
                            else:
                                y = fluxnjy2mag(df_forced[col_y][mask])
                                yerr = flux_err_2_mag_err(
                                    df_forced[col_yerr][mask],
                                    df_forced[col_y][mask].abs(),
                                )

                                if use_absmag and z_obj is not None:
                                    y = absmag(y, distmod=distmod)

                            if use_folded and period is not None:
                                phase = (
                                    np.mod(df_forced[col_mjd][mask], period) / period
                                )
                                for nphase in [0, 1]:
                                    if nphase == 0:
                                        label = "%s (forced)" % band
                                    else:
                                        label = None
                                    ax_i.errorbar(
                                        phase + nphase,
                                        y,
                                        yerr=yerr,
                                        alpha=alpha,
                                        c=colors[band],
                                        linestyle="None",
                                    )
                                    ax_i.scatter(
                                        phase + nphase,
                                        y,
                                        s=sizes_forced[band],
                                        alpha=alpha,
                                        c=colors[band],
                                        marker=markers[band],
                                        label=label,
                                        linestyle="None",
                                    )
                            else:
                                ax_i.errorbar(
                                    df_forced[col_mjd][mask],
                                    y,
                                    yerr=yerr,
                                    alpha=alpha,
                                    c=colors[band],
                                    linestyle="None",
                                )
                                ax_i.scatter(
                                    df_forced[col_mjd][mask],
                                    y,
                                    s=sizes_forced[band],
                                    alpha=alpha,
                                    c=colors[band],
                                    marker=markers[band],
                                    label="%s (forced)" % band,
                                    linestyle="None",
                                )

                if use_folded and ~np.isnan(period):
                    ax_i.set_xlabel("Phase [days]")
                else:
                    ax_i.set_xlabel("MJD")

                if use_flux:
                    if light_type == "diff":
                        ax_i.set_ylabel("Difference flux [nJy]")
                    elif light_type == "sci":
                        ax_i.set_ylabel("Science flux [nJy]")
                    else:
                        ax_i.set_ylabel("")
                    title_aux = d_objs[sid][oid]["lc_kwargs"][
                        "title_exts_flux_" + light_type
                    ]
                    ax_i.set_title(title + title_aux + title_exts)
                else:
                    if light_type == "diff":
                        if use_absmag:
                            ax_i.set_ylabel("Absolute magnitude (from difference flux)")
                        else:
                            ax_i.set_ylabel("Difference magnitude")
                    elif light_type == "sci":
                        if use_absmag:
                            ax_i.set_ylabel("Absolute magnitude")
                        else:
                            ax_i.set_ylabel("Apparent magnitude")
                    else:
                        ax_i.set_ylabel("")
                    title_aux = d_objs[sid][oid]["lc_kwargs"][
                        "title_exts_mag_" + light_type
                    ]
                    ax_i.set_title(
                        title + title_aux + title_exts + " (only positive fluxes shown)"
                    )
                    ax_i.set_ylim(ax_i.get_ylim()[::-1])

                ax_i.legend(loc="center left", bbox_to_anchor=(1, 0.5))

                i += 1

    if mjd_lims is not None:
        for j in fig.axes:
            j.set_xlim([mjd_min, mjd_max])
    else:
        xlims = []
        for j in fig.axes:
            xlims = np.append(xlims, j.get_xlim())
        xlim1 = min(xlims)
        xlim2 = max(xlims)
        for j in fig.axes:
            j.set_xlim([xlim1, xlim2])

    if y_lims is not None:
        for j in fig.axes:
            j.set_ylim(y_lims)

    if show_yr and not use_folded:
        for j in fig.axes:
            xlims = j.get_xlim()
            ylims = j.get_ylim()
            y_yr = (ylims[1] - ylims[0]) * 0.93 + ylims[0]

            mask = (yr["mjd_jan1"] >= xlims[0] - 365.0) & (yr["mjd_jan1"] <= xlims[1])

            for x_yr, label_yr in zip(yr["mjd_jan1"][mask], yr["label"][mask]):
                j.axvline(
                    x=x_yr, color="black", linestyle=":", linewidth=0.7, alpha=0.5
                )
                j.text(
                    x_yr + 0.5 * 365.0,
                    y_yr,
                    label_yr,
                    alpha=0.4,
                    clip_on=True,
                    horizontalalignment="center",
                    fontsize="small",
                )

                j.set_xlim(xlims)

    fig.tight_layout()

    if namefig is not None:
        fig.savefig(namefig, bbox_inches="tight")

    if fontsize is not None:
        plt.rcParams.update({"font.size": 12})
