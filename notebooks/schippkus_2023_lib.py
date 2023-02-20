def generate_synthetic_wave():
    pass


def plot_results_on_map():
    pass


def convert_latlon_to_cartesian(coords):
    """
    Convert an Nx2 array of lon/lat coordinates to cartesian coordinates in [km]

    :param coords: Array of lon/lat coordinates, shape Nx2.
    :type coords: numpy.ndarray
    :return: Array of relative cartesian coordinates, shape Nx2.
    :rtype: numpy.ndarray
    """

    from numpy import array, mean
    from obspy.signal.util import util_geo_km

    coords_cartesian = []
    for sta in coords:
        coords_cartesian.append(
            util_geo_km(mean(coords[:, 0]), mean(coords[:, 1]), sta[0], sta[1])
        )

    coords_cartesian = array(coords_cartesian)
    return coords_cartesian


def parse_station_coordinates_from_inventory(inv_file="../meta/graefenberg.stationxml"):
    from numpy import array
    from obspy import read_inventory

    inv = read_inventory(inv_file)

    coords_latlon = []
    station_coords_sorting = []
    for sta in inv[0]:
        coords_latlon.append([sta.longitude, sta.latitude])
        station_coords_sorting.append(sta.code)

    return array(coords_latlon), station_coords_sorting


def parse_correlations_from_disk(settings):
    from numpy import array
    from obspy import read

    corrs = []
    corr_dir = "../correlations/"
    for sta in settings["station_coords_sorting"]:
        st = read(
            f"{corr_dir}/{settings['master_station'].split('.')[1]}.{sta}.sac",
            format="sac",
        )
        corrs.append(st[0].data)
    return array(corrs)


def select_settings_from_file(label):
    from tomllib import load

    with open("settings.toml", "rb") as f:
        settings = load(f)

    return settings[label]


def compute_source_time_functions(
    settings, sources, times, freqs_excited, apply_excitation_pattern=False
):
    import numpy as np
    from scipy.signal import fftconvolve, ricker
    from tqdm import tqdm

    # print(times.shape, freqs_excited.shape)
    if settings["source_type"] == "microseism":
        all_phases = np.random.uniform(
            0, 2 * np.pi, size=(len(sources), *freqs_excited.shape)
        )

        if settings["phase_mode"] == "constant":
            all_phases = np.repeat(
                [np.random.uniform(0, 2 * np.pi, size=(freqs_excited.shape))],
                repeats=len(sources),
                axis=0,
            )
        # stf = sum over cosine functions
        # one cosine for each excited frequency with a random phase
        # time-frequency combinations don't change
        ft = 2 * np.pi * freqs_excited[:, np.newaxis] * times
        stfs_for_sources = []
        for source_phases in tqdm(all_phases, desc="source terms"):
            # add source's phases and sum over cosines
            stf = np.sum(np.cos(ft + source_phases[:, np.newaxis]), axis=0)
            # all_cosines = [np.cos(ft[idx, :] + p) for idx, p in enumerate(phs)]
            stfs_for_sources.append(stf)

    elif settings["source_type"] == "ricker":
        # excitation_pattern = np.zeros(len(times))
        # excitation_pattern[int(interval*freq):50*int(interval*freq):int(interval*freq)] = 1
        # excitation_pattern[int(interval*freq) + :] = 0
        wavelet = ricker(len(times), 2 * settings["freq"])
        # stfs = [fftconvolve(wavelet, excitation_pattern) for _ in phases]

        if apply_excitation_pattern:
            excitation_pattern = np.zeros(len(times))
            excitation_pattern[
                : settings["interval_count"]
                * int(settings["interval"] * settings["freq"]) : int(
                    settings["interval"] * settings["freq"]
                )
            ] = 1

            stfs_for_sources = [
                fftconvolve(wavelet, excitation_pattern) for _ in sources
            ]

        else:
            stfs_for_sources = [wavelet for _ in sources]

    return np.array(stfs_for_sources)


def compute_synthetic_seismograms_for_sources(station, sources, stfs, freqs, settings):
    import numpy as np
    from obspy.geodetics import gps2dist_azimuth
    from scipy.signal import fftconvolve

    waveforms = []
    for source, stf in zip(sources, stfs):
        dist = gps2dist_azimuth(*station[::-1], *source[::-1])[0] / 1000
        traveltime = dist / settings["medium_velocity"]
        medium_response = np.fft.ifft(np.exp(-1j * 2 * np.pi * freqs * traveltime)).real
        waveform = fftconvolve(medium_response, stf, mode="same")
        # waveform = fftconvolve(waveform, excitation_pattern, mode="same")
        waveform /= np.max(np.abs(waveform))
        waveforms.append(waveform)
    return np.array(waveforms)


def define_isolated_sources(settings):
    import numpy as np

    isolated_sources_per_region = []
    for center_lon, center_lat in zip(settings["center_lons"], settings["center_lats"]):
        isolated_sources = np.array(
            np.random.uniform(low=-1, high=1, size=(settings["n_isolated_sources"], 2))
        )
        isolated_sources[:, 0] += center_lon
        isolated_sources[:, 1] += center_lat
        isolated_sources_per_region.append(isolated_sources)
    return isolated_sources_per_region


def define_boundary_sources(settings):
    import numpy as np
    from pygc import great_circle

    mean_point = np.mean(
        [settings["array_coords_mean"], settings["master_station_coords"]], axis=0
    )

    azs = np.degrees(np.linspace(0, 2 * np.pi, settings["n_boundary_sources"]))[:-1]
    boundary_sources = []
    for az in azs:
        p = great_circle(
            distance=settings["radius"] * 1000,
            azimuth=az,
            latitude=mean_point[1],
            longitude=mean_point[0],
        )
        boundary_sources.append([p["longitude"], p["latitude"]])
    return np.array(boundary_sources)


def compute_synthetic_waveforms():
    pass


def compute_synthetic_correlations(
    boundary_sources, isolated_sources_per_region, settings
):
    import numpy as np
    from scipy.signal import fftconvolve
    from tqdm import tqdm

    times = np.arange(
        0, settings["timelength"] + 1 / settings["freq"], 1 / settings["freq"]
    )
    freqs = np.fft.fftfreq(len(times), 1 / settings["freq"])

    freqs_excited = freqs[
        np.where((freqs > settings["fmin"]) & (freqs < settings["fmax"]))
    ]

    # freqs_excited = np.linspace(settings["fmin"], settings["fmax"], 100)

    # compute source time functions
    stfs_boundary_sources = compute_source_time_functions(
        settings,
        boundary_sources,
        times,
        freqs_excited,
        apply_excitation_pattern=settings["apply_excitation_pattern_for_boundary"],
    )

    stfs_isolated_sources_per_region = []
    for isolated_sources in isolated_sources_per_region:
        stfs_isolated_sources = compute_source_time_functions(
            settings,
            isolated_sources,
            times,
            freqs_excited,
            apply_excitation_pattern=settings["apply_excitation_pattern_for_isolated"],
        )
        stfs_isolated_sources_per_region.append(stfs_isolated_sources)

    # isolated_source_mean = np.mean(isolated_sources, axis=0)
    # compute waveforms for array stations
    boundary_source_waveforms_per_array_station = []
    for station in tqdm(settings["station_coords"], desc="waveforms"):
        boundary_source_waveforms_per_array_station.append(
            compute_synthetic_seismograms_for_sources(
                station=station,
                sources=boundary_sources,
                stfs=stfs_boundary_sources,
                freqs=freqs,
                settings=settings,
            )
        )

    isolated_source_waveforms_per_array_station_per_region = []
    for isolated_sources, stfs_isolated_sources in zip(
        isolated_sources_per_region, stfs_isolated_sources_per_region
    ):
        isolated_source_waveforms_per_array_station = []
        for station in tqdm(settings["station_coords"], desc="waveforms"):
            isolated_source_waveforms_per_array_station.append(
                np.mean(
                    compute_synthetic_seismograms_for_sources(
                        station=station,
                        sources=isolated_sources,
                        stfs=stfs_isolated_sources,
                        freqs=freqs,
                        settings=settings,
                    ),
                    axis=0,
                )
            )
        isolated_source_waveforms_per_array_station_per_region.append(
            isolated_source_waveforms_per_array_station
        )

    # compute waveforms for master station
    boundary_source_waveforms_master_station = compute_synthetic_seismograms_for_sources(
        station=settings["master_station_coords"],
        sources=boundary_sources,
        stfs=stfs_boundary_sources,
        freqs=freqs,
        settings=settings,
    )

    isolated_source_waveforms_master_station_per_region = []
    for isolated_sources, stfs_isolated_sources in zip(
        isolated_sources_per_region, stfs_isolated_sources_per_region
    ):
        isolated_source_waveforms_master_station = compute_synthetic_seismograms_for_sources(
            station=settings["master_station_coords"],
            sources=isolated_sources,
            stfs=stfs_isolated_sources,
            freqs=freqs,
            settings=settings,
        )
        isolated_source_waveforms_master_station_per_region.append(
            isolated_source_waveforms_master_station
        )
    isolated_source_waveforms_master_station_per_region = np.mean(
        isolated_source_waveforms_master_station_per_region, axis=0
    )

    mean_isolated_source_waveforms_master_station = np.mean(
        isolated_source_waveforms_master_station, axis=0
    )

    # correlate boundary source waveforms
    boundary_source_correlations = []
    for waveform in tqdm(
        boundary_source_waveforms_per_array_station,
        desc="correlations boundary sources",
    ):
        correlations = fftconvolve(
            waveform,
            boundary_source_waveforms_master_station[:, ::-1],
            mode="same",
            axes=1,
        )
        correlations = np.mean(correlations, axis=0)
        boundary_source_correlations.append(correlations)
    boundary_source_correlations = np.array(boundary_source_correlations)

    # correlate isolated sources waveforms
    isolated_source_correlations_per_region = []
    for (isolated_source_waveforms_per_array_station, region_weight) in zip(
        isolated_source_waveforms_per_array_station_per_region,
        settings["region_weights"],
    ):
        isolated_source_correlations = []
        for waveform in tqdm(
            isolated_source_waveforms_per_array_station,
            desc="correlations isolated sources",
        ):
            # print(waveform.shape)
            # print(mean_isolated_source_waveforms_master_station.shape)
            correlations = fftconvolve(
                waveform,
                mean_isolated_source_waveforms_master_station[::-1],
                # isolated_source_waveforms_master_station[:, ::-1],
                mode="same",
                # axes=1,
            )
            # correlations = np.mean(correlations, axis=0)
            isolated_source_correlations.append(correlations)
        isolated_source_correlations = np.array(isolated_source_correlations)
        isolated_source_correlations /= np.max(np.abs(isolated_source_correlations))
        isolated_source_correlations *= region_weight
        isolated_source_correlations_per_region.append(isolated_source_correlations)

    isolated_source_correlations = np.mean(
        isolated_source_correlations_per_region, axis=0
    )

    # isolated_source_correlations /= np.max(np.abs(isolated_source_correlations))
    boundary_source_correlations /= np.max(np.abs(boundary_source_correlations))

    # -- amplitude scaling
    isolated_source_correlations *= settings["contribution_amplitude_ratio"]

    corrs = boundary_source_correlations + isolated_source_correlations
    # corrs = boundary_source_correlations

    return corrs


def beamforming(data, slowness_space, settings, lapsetimes, windows, coordinates):

    import numpy as np

    # estimate best starttime and endtime points
    starttime_idxs = [
        np.argmin(np.abs(lapsetimes - starttime)) for starttime in windows
    ]
    endtime_idxs = [
        np.argmin(np.abs(lapsetimes - (starttime + settings["window_length"])))
        for starttime in windows
    ]

    # compute expected frequencies, and limit to range of interest
    f = np.fft.fftfreq(
        settings["window_length"] * settings["freq"], d=1 / settings["freq"]
    )
    omega = 2 * np.pi * f
    omega_limited = omega[np.where((f > settings["fmin"]) & (f < settings["fmax"]))]

    # cut timewindows from correlations
    cut_corrs_all = np.zeros(
        (len(windows), data.shape[0], settings["window_length"] * settings["freq"])
    )
    for window_idx, (starttime_idx, endtime_idx) in enumerate(
        zip(starttime_idxs, endtime_idxs)
    ):
        cut_corrs_all[window_idx] = data[:, starttime_idx:endtime_idx]

    # compute spectra and limit to frequency band of interest
    data_spectra_all = np.fft.fft(cut_corrs_all, axis=2)
    data_spectra_all_limited = data_spectra_all[
        :, :, np.where((f > settings["fmin"]) & (f < settings["fmax"]))[0]
    ]

    zero_spectra = np.zeros(omega_limited.shape).astype(complex)

    # compute all traveltimes for all slowness - station combinations
    # relative to the reference point, which is arbitrary but must be reasonably nearby
    #
    # einsum dimension labels
    # - n: station dimension, number of stations
    # - s: slowness dimension, number of slowness grid points
    # - x: spatial dimensions, here 2D (only surface coordinates)
    traveltimes_space = np.einsum(
        "nx, sx -> sn", settings["reference_point"] - coordinates, slowness_space
    )

    # compute synthetic spectra, i.e., plane waves, for traveltimes above
    # other names: replica vectors, Green's functions
    #
    # einsum dimension labels
    # - w: frequency dimension, number of frequencies in band [fmin, fmax]
    # - s: slowness dimension, number of slowness grid points
    # - n: station dimension, number of stations
    synth_spectra_space = np.exp(
        -1j * np.einsum("w, sn -> snw", omega_limited, traveltimes_space)
    )

    # cross-spectral density matrix K for each beamforming window
    # contains cross correlations of all input signals
    #
    # einsum dimension labels
    # - t: timewindow dimension, number of beampower windows
    # - i: first station dimension, number of stations
    # - j: second station dimension, number of stations
    # - w: frequency dimension, number of frequencies in band [fmin, fmax]
    K = np.einsum(
        "tiw, tjw -> tijw", data_spectra_all_limited, np.conj(data_spectra_all_limited)
    )
    # K = np.einsum("iw,jw->ijw", data_spectra, np.conj(data_spectra), optimize=True)

    # remove auto-correlations, i.e., energy-scaling of beampowers
    # this yields beampowers in range [-max, max]
    for idx in range(K.shape[1]):
        K[:, idx, idx, :] = zero_spectra

    # conventional beamformer
    # B = Σ_ω Σ_i Σ_j s^*_i(ω) K_ij(ω) s_j(ω)
    # for more details, see e.g., Schippkus et al. 2022

    # einsum dimension labels
    # - t: timewindow dimension, number of beampower windows
    # - s: slowness dimension, number of slowness grid points
    # - i: first station dimension, number of stations
    # - j: second station dimension, number of stations
    # - w: frequency dimension, number of frequencies in band [fmin, fmax]

    beampowers_for_time = np.einsum(
        "siw, tijw, sjw -> ts", np.conj(synth_spectra_space), K, synth_spectra_space
    ).real
    return np.array(beampowers_for_time)


def extract_beampower_peaks(beampowers_for_time, slows, azs):
    import numpy as np

    peaks = []
    for bp in beampowers_for_time:
        bp = bp.reshape(len(azs), len(slows))
        peak_idxs = np.where(bp == np.max(bp))
        vel = 1 / slows[peak_idxs[1][0]]
        baz = 90 - np.degrees(azs[peak_idxs[0][0]])
        amp = np.max(bp)
        peaks.append([vel, baz, amp])
    return np.array(peaks)


def compute_correlation_with_beams(
    windows, peaks, lapsetimes, settings, corrs_filt, coords_cartesian
):
    import numpy as np

    # define frequencies for beamforming windows
    freqs = np.fft.fftfreq(
        settings["window_length"] * settings["freq"], 1 / settings["freq"]
    )
    omega = 2 * np.pi * freqs

    corr_coeffs_for_windows = []
    for starttime, (vel, baz, amp) in zip(windows, peaks):

        # extract time window from correlations
        starttime_idx = np.argmin(np.abs(lapsetimes - starttime))
        endtime_idx = np.argmin(
            np.abs(lapsetimes - (starttime + settings["window_length"]))
        )

        cut_corrs = corrs_filt[:, starttime_idx:endtime_idx]
        cut_corrs_spectra = np.fft.fft(cut_corrs, axis=1)

        # compute relative traveltimes for best-beam slowness
        slowness = (
            np.cos(np.radians(-baz + 90)) * 1 / vel,
            np.sin(np.radians(-baz + 90)) * 1 / vel,
        )

        traveltimes = np.dot(coords_cartesian - settings["reference_point"], slowness)

        # plane waves
        replica_vectors = np.exp(-1j * np.einsum("w, n->nw", omega, traveltimes))

        # apply time shifts corresponding to best-beam plane-wave
        # i.e., project onto replica vectors
        corrs_shift = np.fft.ifft(cut_corrs_spectra * replica_vectors).real

        beam = np.mean(corrs_shift, axis=0)
        corr_coeffs = []
        for c in corrs_shift:
            corr_coeffs.append(np.corrcoef(c, beam)[0, 1])

        corr_coeffs_for_windows.append(corr_coeffs)
    return np.array(corr_coeffs_for_windows)
