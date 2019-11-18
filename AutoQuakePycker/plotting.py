#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Routines for plotting the results of AutoQuakePycker
"""
import matplotlib.pyplot as plt
import numpy as np


def plot_sta_results(sta, st_sta, runtype, orig_time, pdf, fig, cmps, cfg):
    """
    Plot picking results for each station

    Parameters
    ----------
    sta : Dictionary of station pick results
    st_sta : ObsPy stream containing waveform data for station
    runtype :
    tr : ObsPy trace object

    Returns
    ----------
    ev_new : ObsPy event object containing relocated event
    """
    co = ["blue", "orange"]
    ct = ["green", "red"]
    for n_tr, tr in enumerate(st_sta):
        cmp = tr.stats.channel[2:3]
        if cmp in cmps["P"]:
            phase = "P"
            nphs = 0
        elif cmp in cmps["S"]:
            phase = "S"
            nphs = 1
        sta_ = sta[phase]
        try:
            sta_[cmp]
        except:
            continue
        ax1 = plt.subplot(4, len(st_sta), n_tr+1)
        ax2 = plt.subplot(4, len(st_sta), n_tr+len(st_sta)+1, sharex=ax1)
        ax3 = plt.subplot(4, len(st_sta), n_tr+2*(len(st_sta))+1, sharex=ax1)
        ax4 = plt.subplot(4, len(st_sta), n_tr+3*(len(st_sta))+1, sharex=ax1)
        if len(sta_[cmp]["times"]) != len(sta_[cmp]["f1_mean_smooth"]):
            continue
        ax2.plot(sta_[cmp]["times"], sta_[cmp]["f1_mean_smooth"],
                 label="F$_4$")
        for n_filt, filt in enumerate(sta_[cmp]["filtwins_check"]):
            ax1 = ax1.twinx()
            ax1.plot(sta_[cmp]["times"], sta_[cmp]["tr_results"][n_filt],
                     linewidth=0.2, color="black")
            ax1.set_yticklabels([])
            ax1.set_ylim([
                    -1*np.max(np.abs(sta_[cmp]["tr_results"][n_filt])),
                    np.max(np.abs(sta_[cmp]["tr_results"][n_filt]))
            ])
            for n_kurt, kurt in enumerate(cfg.picking.KURT_WINS):
                ax3.plot(sta_[cmp]["times"],
                         sta_[cmp]["f3_results"][n_filt][n_kurt],
                         linewidth=0.7, linestyle="--")
        p = sta["best_obs"][0]
        for n in range(0, len(sta["best_obs"])):
            p = sta["best_obs"][n]
            if n == 0:
                phase = "P"
            else:
                phase = "S"
            if p:
                if p[2] == tr.stats.channel:
                    ax1.text(p[0]-orig_time, ax1.get_ylim()[1],
                             "{:1g}".format(p[3]),
                             horizontalalignment="center", bbox=dict(
                                     facecolor='none', edgecolor='red'))
                if cmp in cmps[phase]:
                    ax1.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10,
                                label="Best Kurtosis {:}-pick".format(phase))
                    ax2.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax3.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax4.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax1.set_xlim([
                            p[0]-orig_time-cfg.plotting.PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+cfg.plotting.PLOT_WIND_PICK_DIFF
                                  ])
                    ax2.set_xlim([
                            p[0]-orig_time-cfg.plotting.PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+cfg.plotting.PLOT_WIND_PICK_DIFF
                                  ])
                    ax3.set_xlim([
                            p[0]-orig_time-cfg.plotting.PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+cfg.plotting.PLOT_WIND_PICK_DIFF
                                  ])
                    ax4.set_xlim([
                            p[0]-orig_time-cfg.plotting.PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+cfg.plotting.PLOT_WIND_PICK_DIFF
                                  ])

        for phase in ["P", "S"]:
            if cmp in cmps[phase]:
                ax1.axvline(sta["pred"][phase]-orig_time,
                            label="Theor. {:}-arrival".format(phase),
                            color=ct[nphs], linestyle="--")
        for n_sw, sw in enumerate(cfg.picking.CF_MEAN_SMOOTH_WIND):
            ax4.plot(sta_[cmp]["times"], sta_[cmp]["f4_all"][n_sw],
                     label="F$_4$ {:}".format(sw))
        ax4_ = ax3.twinx()
        ax3.plot(sta_[cmp]["times"], sta_[cmp]["f3_mean_smooth"][0],
                 color="red", linewidth=2, label="Smoothed mean CF")
        ax1.text(0.04, 0.96, "Waveforms {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax1.transAxes,
                 bbox=dict(facecolor="white", edgecolor="black",
                           boxstyle="round", alpha=0.5))
        ax2.text(0.04, 0.96, "Raw mean Kurtosis {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax2.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        ax3.text(0.04, 0.96, "Kurtosis CF {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax3.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        ax4.text(0.04, 0.96, "F$_4$ transform {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax4.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        plt.suptitle("{:}.{:}".format(
                st_sta[n_tr].stats.network, st_sta[n_tr].stats.station),
                fontsize=16)
        ax1.grid()
        ax1.legend(loc="lower left")
        ax2.grid()
        ax2.legend(loc="lower left")
        ax3.grid()
        ax3.legend()
        ax4.grid()
        h1, l1 = ax4.get_legend_handles_labels()
        h2, l2 = ax4_.get_legend_handles_labels()
        ax4.legend(h1+h2, l1+l2, loc="upper right")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.subplots_adjust(wspace=0.1)
    pdf.savefig(fig)
    # plt.show()
    plt.cla()
    plt.clf()


def plot_record_section(st_sta):
    """
    TODO
    """
    fig = plt.figure()
    st_sta.plot(type="section", orientation="horizontal",
                vertical_scaling_range=5e3, fig=fig)
#    for tr in st_sta:
#    plt.scatter([30], [60])
