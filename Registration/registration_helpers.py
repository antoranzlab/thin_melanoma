"""Coordinate replay and optional fitting for the supplied VALIS registrations."""
import importlib.metadata
import json
import os
import random

import numpy as np
import pandas as pd


def start_valis(seed=2026):
    os.environ.setdefault('CUDA_VISIBLE_DEVICES', '-1')
    if importlib.metadata.version('valis-wsi') != '1.1.0':
        raise RuntimeError('Use VALIS 1.1.0, the version used for these registrations.')
    random.seed(seed)
    np.random.seed(seed)
    import cv2
    cv2.setRNGSeed(seed)
    from valis import registration
    registration.init_jvm()
    return registration


def get_registrar(config, data_dir, output_dir, refit=False):
    registration = start_valis()
    if not refit:
        return registration.load_registrar(str(data_dir / config['registrar']))
    from valis import feature_detectors
    image_dir = data_dir / config['images']
    if not (image_dir / 'xenium.tif').is_file():
        raise FileNotFoundError(image_dir / 'xenium.tif')
    fit_dir = output_dir / 'fit'
    fit_dir.mkdir(parents=True, exist_ok=False)
    settings = config['fit'].copy()
    detector = settings.pop('feature_detector_cls', None)
    if detector:
        settings['feature_detector_cls'] = getattr(feature_detectors, detector)
    registrar = registration.Valis(src_dir=str(image_dir), dst_dir=str(fit_dir), reference_img_f=str(image_dir / 'xenium.tif'), **settings)
    registrar.register(**config['register'])
    return registrar


def transform_coordinates(frame, modality, config, registrar=None):
    frame = frame.copy()
    if modality == 'xenium':
        frame['sample_id'] = frame['sample_id'].astype(str)
        frame = frame.loc[frame['sample_id'] == config['xenium_acquisition']].copy()
        offset = config['xenium_crop_x_um']
        if offset:
            frame = frame.loc[frame['x'] > offset].copy()
            frame['x'] = frame['x'] - offset
        points = frame[['x', 'y']].to_numpy() * (0.2125 / 0.6525)
        points = points * (1 / 0.2125)
    elif modality == 'milan':
        for column in ['X', 'Y']:
            frame[column] = frame[column].astype(str).str.replace(',', '.', regex=False).astype(float)
        points = registrar.get_slide('milan').warp_xy(frame[['X', 'Y']].to_numpy() * 1, non_rigid=True)
        frame['sample_id'] = frame['slide_id'] + frame['scan_region']
    elif modality == 'visium':
        frame = frame.loc[frame['in_tissue'] == 1].copy()
        points = frame[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy().copy()
        points = (points - np.asarray(config['visium_crop_xy'])) * (0.18 / 0.6525)
        points = registrar.get_slide('visium').warp_xy(points, non_rigid=True)
        frame['sample_id'] = config['visium_source_sample']
    else:
        raise ValueError(modality)
    if frame.empty or not np.isfinite(points).all():
        raise ValueError(f'{modality}: empty selection or non-finite coordinates')
    frame['x_registered'] = points[:, 0]
    frame['y_registered'] = points[:, 1]
    frame['source_sample_id'] = frame['sample_id']
    frame['sample_id'] = config['patient']
    frame['patient'] = config['patient']
    return frame


def read_coordinates(config, data_dir, modality):
    sep = ';' if modality == 'milan' else ','
    return pd.read_csv(data_dir / config['inputs'][modality], sep=sep)


def compare_coordinates(actual, baseline, modality, tolerance_px=1e-6):
    keys = {'milan': ['OID', 'slide_id', 'scan_region'], 'xenium': ['cell'], 'visium': ['barcode']}[modality]
    columns = ['x_registered', 'y_registered']
    actual = actual.copy()
    baseline = baseline.copy()
    for key in keys:
        actual[key] = actual[key].astype(str)
        baseline[key] = baseline[key].astype(str)
    if actual.duplicated(keys).any() or baseline.duplicated(keys).any():
        raise ValueError(f'{modality}: duplicate IDs')
    a = actual.set_index(keys)[columns].sort_index()
    b = baseline.set_index(keys)[columns].sort_index()
    if not a.index.equals(b.index):
        raise ValueError(f'{modality}: baseline and replay have different IDs')
    distances = np.linalg.norm(a.to_numpy() - b.to_numpy(), axis=1)
    if not len(distances) or not np.isfinite(distances).all():
        raise ValueError(f'{modality}: empty or non-finite comparison')
    maximum = float(distances.max())
    return {'modality': modality, 'rows': len(a), 'median_px': float(np.median(distances)), 'max_px': maximum, 'tolerance_px': tolerance_px, 'passed': maximum <= tolerance_px}


def check_and_save(frames, config, data_dir, output_dir):
    reports = []
    for modality, frame in frames.items():
        baseline_path = data_dir / config['baselines'][modality]
        baseline = pd.read_csv(baseline_path)
        report = compare_coordinates(frame, baseline, modality)
        report['baseline'] = str(baseline_path.resolve())
        reports.append(report)
    output_dir.mkdir(parents=True, exist_ok=False)
    pd.DataFrame(reports).to_csv(output_dir / 'coordinate_comparison.csv', index=False)
    if not all(report['passed'] for report in reports):
        raise ValueError('Coordinates differ from the selected baseline. See coordinate_comparison.csv; no coordinate files were saved.')
    for modality, frame in frames.items():
        frame.to_csv(output_dir / f'{modality}_{config["patient"]}.csv', index=False)
    (output_dir / 'run_config.json').write_text(json.dumps(config, indent=2) + '\n')
