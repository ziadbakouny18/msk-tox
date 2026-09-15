#!/usr/bin/env python3
"""
Generate steroid ratio cache — v4.

FILE STRUCTURE (verified against the real CSVs)
-----------------------------------------------
Every llm84k_<tox>_grade<k>_<date>.csv shares an IDENTICAL LOT skeleton:
the same 127,417 rows, same order, same mrn / lot / lot_start / lot_end /
lot_end_pad / contains_* columns. The files differ in exactly four columns:
ae_date, ae_toxicity, ae_grade, t_ae.

The AE overlay holds the FIRST event meeting that file's threshold, not the
patient's max grade. A patient who went grade 1 -> grade 2 appears as grade 1
(earlier date) in the grade0 file and grade 2 (later date) in the grade2 file.
So the threshold files are NOT redundant: reading grade0 alone under-grades
every patient who escalated.

Therefore:
  * The skeleton is read ONCE (from any file) and asserted identical elsewhere.
  * AE columns are read per (toxicity, threshold) and stacked as overlays.
  * grade    = max(ae_grade) across all overlays for that (mrn, toxicity)
  * ref_date = earliest ae_date at that max grade
  * Grades 1 and 4 both survive (95 grade-4 pneumonitis patients exist; a
    threshold-ladder derivation would fold them into grade 3).

OTHER DATA FACTS BAKED IN
-------------------------
  * lot_end_pad == lot_end + 180 days, exactly, on every row. 32% of AEs fall
    after lot_end but 100% fall within lot_end_pad, so ICI attribution uses
    lot_end_pad and there is no separate lookback constant.
  * ae_grade is NaN (not 0) for non-event rows; ~0.7% of rows are events.
  * t_ae is +inf on every non-event row. Never aggregate it.
  * contains_immuno is a real numpy bool in these files, but to_bool() is kept
    as a guard in case a regenerated file lands as object strings.

All steroid-calculation corrections from the review are carried over, tagged
[FIX-n]. Outputs: steroid_ratio_cache_LOT.pkl
"""

import re
import pickle
import warnings
from pathlib import Path
from datetime import timedelta

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)

# ============================================================================
# CONFIG
# ============================================================================

ROOT = Path(__file__).parent.parent.resolve()
FIGURES = ROOT.parent.parent
DATA = FIGURES / "figures_data" / "figure 1" / "data"
LOT_DIR = DATA / "OneDrive_1_8-7-2026"

TOXICITIES = [
    'adrenal_insufficiency', 'colitis', 'pneumonitis',
    'liver_toxicity', 'hypothyroidism', 'hyperthyroidism',
]
WINDOW_DAYS = 60

# ICI is defined identically for every grade: ref_date falls inside an
# ICI-containing LOT's [lot_start, lot_end_pad] window. The 180-day pad already
# encodes the attribution window, so no extra lookback is applied.
LOT_END_COL = 'lot_end_pad'      # 'lot_end_pad' | 'lot_end'

# Grade-0 controls have no AE date. Anchor them to the midpoint of an
# ICI-containing LOT so they sit relative to ICI the way cases do.
GRADE0_ANCHOR = 'ici_lot'          # 'ici_lot' | 'full_timeline'

EPSILON = 0.01

INCLUDE_FLUDROCORTISONE = True
MAX_COURSE_DAYS = 365
DEFAULT_DOSE_MG = 40.0

# Prednisone-equivalent multipliers: raw mg * factor = prednisone-equivalent mg.
# Anchor: hydrocortisone 20 = prednisone 5 = methylprednisolone 4 = dex 0.75.
STEROID_EQUIVALENTS = {
    'cortisone':          5.0 / 25.0,   # 0.20
    'hydrocortisone':     5.0 / 20.0,   # 0.25
    'prednisone':         1.0,
    'prednisolone':       1.0,
    'triamcinolone':      5.0 / 4.0,    # 1.25
    'methylprednisolone': 5.0 / 4.0,    # 1.25
    'betamethasone':      5.0 / 0.6,    # 8.33
    'dexamethasone':      5.0 / 0.75,   # 6.67
    'fludrocortisone':    5.0 / 2.0,    # 2.50 (GC potency; correct as written)
}
if not INCLUDE_FLUDROCORTISONE:
    STEROID_EQUIVALENTS.pop('fludrocortisone')

# [FIX-1] Longest name first.
STEROID_KEYS = sorted(STEROID_EQUIVALENTS, key=len, reverse=True)

NON_SYSTEMIC_KEYWORDS = [
    'topical', 'cream', 'ointment', 'lotion', 'gel', 'ophthalmic', 'otic',
    'eye', 'ear drop', 'inhaler', 'inhalation', 'nebulizer', 'nasal',
    'rectal', 'suppository', 'enema', 'patch', 'intra-articular',
    'intraarticular',
]

RX_FREQ_COL = None
RX_FREQ_CANDIDATES = ['FREQUENCY', 'FREQ', 'FREQ_DESC', 'ADMIN_FREQ', 'SIG',
                      'DOSE_FREQ', 'FREQUENCY_DESC']

THRESHOLD_RE = re.compile(r'grade(\d+)', re.IGNORECASE)

DROPS = {}


def record_drop(step, before, after):
    lost = before - after
    DROPS[step] = lost
    pct = (100.0 * lost / before) if before else 0.0
    print(f"    {step}: {before:,} -> {after:,}  (dropped {lost:,}, {pct:.1f}%)")


# ============================================================================
# HELPERS
# ============================================================================

def standardize_mrn(mrn):
    if pd.isna(mrn):
        return None
    digits = ''.join(filter(str.isdigit, str(mrn).strip()))
    return digits.zfill(8) if digits else None


def to_bool(series):
    """[FIX-9] `.astype(bool)` on object strings makes 'False' -> True."""
    if series.dtype == bool:
        return series
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float) > 0
    s = series.astype(str).str.strip().str.lower()
    true_vals = {'true', 't', 'yes', 'y', '1', '1.0'}
    false_vals = {'false', 'f', 'no', 'n', '0', '0.0', 'nan', 'none', '', '(null)'}
    unknown = set(s.unique()) - true_vals - false_vals
    if unknown:
        raise ValueError(f"Unrecognised boolean values: {sorted(unknown)[:10]}")
    return s.isin(true_vals)


def identify_steroid(med_name):
    if pd.isna(med_name):
        return None
    m = str(med_name).lower()
    for excl in NON_SYSTEMIC_KEYWORDS:
        if excl in m:
            return None
    for s in STEROID_KEYS:                 # [FIX-1] longest first
        if s in m:
            return s
    return None


_FREQ_PATTERNS = [
    (r'\bq\s*1\s*h\b', 24.0), (r'\bq\s*2\s*h\b', 12.0), (r'\bq\s*3\s*h\b', 8.0),
    (r'\bq\s*4\s*h\b', 6.0),  (r'\bq\s*6\s*h\b', 4.0),  (r'\bq\s*8\s*h\b', 3.0),
    (r'\bq\s*12\s*h\b', 2.0), (r'\bq\s*24\s*h\b', 1.0), (r'\bq\s*48\s*h\b', 0.5),
    (r'\bqid\b', 4.0), (r'\btid\b', 3.0), (r'\bbid\b', 2.0),
    (r'\bqd\b', 1.0), (r'\bdaily\b', 1.0), (r'\bevery day\b', 1.0),
    (r'\bfour times (a |per )?day\b', 4.0),
    (r'\bthree times (a |per )?day\b', 3.0),
    (r'\btwice (a |per )?day\b', 2.0),
    (r'\bonce (a |per )?day\b', 1.0),
    (r'\bqod\b', 0.5), (r'\bevery other day\b', 0.5),
    (r'\bweekly\b', 1.0 / 7.0),
]


def parse_doses_per_day(freq_text):
    """[FIX-4] v1 computed dose * days, i.e. assumed once-daily everywhere."""
    if pd.isna(freq_text):
        return 1.0
    t = str(freq_text).lower()
    for pattern, per_day in _FREQ_PATTERNS:
        if re.search(pattern, t):
            return per_day
    return 1.0


def expand_to_days(df):
    """[FIX-2] Explode half-open [start, end_ex) into one row per calendar day.

    Day count is (end_ex - start).days with NO +1. Dates are normalized to
    midnight upstream [FIX-11], so this is exact calendar-day arithmetic.
    """
    df = df[df['end_ex'] > df['start']].copy()
    if len(df) == 0:
        return pd.DataFrame(columns=['MRN', 'steroid', 'day', 'mg_pred', 'source'])

    n_days = (df['end_ex'] - df['start']).dt.days.clip(upper=MAX_COURSE_DAYS).values
    total = int(n_days.sum())
    if total == 0:
        return pd.DataFrame(columns=['MRN', 'steroid', 'day', 'mg_pred', 'source'])

    row_idx = np.repeat(np.arange(len(df)), n_days)
    starts = np.cumsum(n_days) - n_days
    offsets = np.arange(total) - np.repeat(starts, n_days)

    return pd.DataFrame({
        'MRN':     df['MRN'].values[row_idx],
        'steroid': df['steroid'].values[row_idx],
        'day':     df['start'].values[row_idx] + offsets.astype('timedelta64[D]'),
        'mg_pred': df['mg_pred_per_day'].values[row_idx],
        'source':  df['source'].values[row_idx],
    })


# ============================================================================
# 1. LOAD SKELETON ONCE + AE OVERLAYS PER (TOXICITY, THRESHOLD)
# ============================================================================

print("=" * 74)
print("1. LOADING LOT SKELETON AND AE OVERLAYS")
print("=" * 74)

SKEL_COLS = ['mrn', 'lot', 'lot_start', 'lot_end', 'lot_end_pad', 'contains_immuno']
AE_COLS = ['ae_date', 'ae_grade', 'ae_toxicity']

file_index = []
for tox in TOXICITIES:
    files = sorted(LOT_DIR.glob(f"llm84k_{tox}_grade*.csv"))
    if not files:
        raise FileNotFoundError(f"No LOT files for '{tox}' in {LOT_DIR}")
    for f in files:
        m = THRESHOLD_RE.search(f.name)
        if not m:
            raise ValueError(f"Cannot parse threshold from filename: {f.name}")
        file_index.append((tox, int(m.group(1)), f))
        print(f"  {f.name:56s} tox={tox:22s} threshold>={m.group(1)}")

# --- Skeleton, read once ---------------------------------------------------
skeleton = pd.read_csv(file_index[0][2], usecols=SKEL_COLS, low_memory=False)
print(f"\n  Skeleton: {len(skeleton):,} LOT rows from {file_index[0][2].name}")

skeleton['mrn'] = skeleton['mrn'].apply(standardize_mrn)
for col in ['lot_start', 'lot_end', 'lot_end_pad']:
    skeleton[col] = pd.to_datetime(skeleton[col], errors='coerce').dt.normalize()  # [FIX-11]
skeleton['ici_status'] = to_bool(skeleton['contains_immuno'])                      # [FIX-9]
skeleton['lot_end_eff'] = skeleton[LOT_END_COL]

n0 = len(skeleton)
skeleton = skeleton[skeleton['mrn'].notna()].copy()
record_drop("skeleton: unparseable MRN", n0, len(skeleton))
print(f"  Unique patients: {skeleton['mrn'].nunique():,}")

pad_gap = (skeleton['lot_end_pad'] - skeleton['lot_end']).dt.days
print(f"  lot_end_pad - lot_end: min={pad_gap.min()} max={pad_gap.max()} days "
      f"(using '{LOT_END_COL}' for ICI attribution)")

# --- AE overlays -----------------------------------------------------------
# Positional alignment is only valid because the skeletons are byte-identical.
# Assert it rather than trust it: one file regenerated on a different day would
# break this silently.
overlays = []
for tox, thr, f in file_index:
    chk = pd.read_csv(f, usecols=['mrn', 'lot'], low_memory=False)
    if len(chk) != n0 or not (chk['lot'].values == pd.read_csv(
            file_index[0][2], usecols=['lot'], low_memory=False)['lot'].values).all():
        raise AssertionError(
            f"{f.name} does not share the reference LOT skeleton. "
            "Positional overlay is unsafe — merge on (mrn, lot) instead.")

    ae = pd.read_csv(f, usecols=AE_COLS, low_memory=False)
    ae['mrn'] = chk['mrn'].apply(standardize_mrn).values
    ae['lot'] = chk['lot'].values
    ae['toxicity'] = tox
    ae['file_threshold'] = thr
    ae['ae_date'] = pd.to_datetime(ae['ae_date'], errors='coerce').dt.normalize()
    ae['ae_grade'] = pd.to_numeric(ae['ae_grade'], errors='coerce')
    # ae_grade is NaN (not 0) on non-event rows; events are the notna rows.
    overlays.append(ae[ae['ae_date'].notna() & ae['ae_grade'].notna()].copy())

events = pd.concat(overlays, ignore_index=True)
events['ae_grade'] = events['ae_grade'].astype(int)
n0 = len(events)
events = events[events['mrn'].notna()].copy()
record_drop("events: unparseable MRN", n0, len(events))

print(f"\n  Event rows across all overlays: {len(events):,}")
print(events.groupby(['toxicity', 'file_threshold'])['ae_grade']
            .agg(['size', 'min', 'max']).rename(columns={'size': 'rows'}))


# ============================================================================
# 2. PATIENT-LEVEL GRADE AND REFERENCE DATE
# ============================================================================

print("\n" + "=" * 74)
print("2. PATIENT-LEVEL AGGREGATION")
print("=" * 74)

# grade = max across ALL overlays. Reading grade0 alone would report the
# first-event grade and miss every escalation.
max_grade = (events.groupby(['mrn', 'toxicity'])['ae_grade'].max()
                   .rename('grade').reset_index())

# ref_date = EARLIEST ae_date at that max grade. [FIX-5] (not .iloc[0])
ev = events.merge(max_grade, on=['mrn', 'toxicity'])
ev = ev[ev['ae_grade'] == ev['grade']]
ref = (ev.groupby(['mrn', 'toxicity'])['ae_date'].min()
         .rename('ref_date').reset_index())
cases = max_grade.merge(ref, on=['mrn', 'toxicity'])

first_grade = (events.sort_values('ae_date')
                     .drop_duplicates(subset=['mrn', 'toxicity'], keep='first')
                     .set_index(['mrn', 'toxicity'])['ae_grade'])
esc = cases.set_index(['mrn', 'toxicity'])['grade'] > first_grade.reindex(
    cases.set_index(['mrn', 'toxicity']).index)
print(f"  Case patient-toxicities: {len(cases):,}")
print(f"  Regraded upward by the threshold overlays: {int(esc.fillna(False).sum()):,}")
print("\n  Max-grade distribution:")
print(cases.pivot_table(index='toxicity', columns='grade', values='mrn',
                        aggfunc='count', fill_value=0))

# --- Controls: patient-toxicities with no event in ANY overlay -------------
pts = skeleton[['mrn']].drop_duplicates()
all_pt = pts.merge(pd.DataFrame({'toxicity': TOXICITIES}), how='cross')
controls = all_pt.merge(cases[['mrn', 'toxicity']], on=['mrn', 'toxicity'],
                        how='left', indicator=True)
controls = controls[controls['_merge'] == 'left_only'][['mrn', 'toxicity']].copy()
print(f"\n  Control (grade 0) patient-toxicities: {len(controls):,}")

# Anchor controls to the midpoint of their ICI treatment span. [FIX-6]
ici_span = (skeleton[skeleton['ici_status']].groupby('mrn')
            .agg(ici_lo=('lot_start', 'min'), ici_hi=('lot_end', 'max')))
all_span = skeleton.groupby('mrn').agg(all_lo=('lot_start', 'min'),
                                       all_hi=('lot_end', 'max'))
controls = (controls.merge(ici_span, on='mrn', how='left')
                    .merge(all_span, on='mrn', how='left'))
ici_mid = controls['ici_lo'] + (controls['ici_hi'] - controls['ici_lo']) / 2
all_mid = controls['all_lo'] + (controls['all_hi'] - controls['all_lo']) / 2
if GRADE0_ANCHOR == 'ici_lot':
    controls['ref_date'] = ici_mid.fillna(all_mid)
    print(f"  Controls anchored to an ICI LOT: {ici_mid.notna().sum():,}"
          f" / {len(controls):,}")
else:
    controls['ref_date'] = all_mid
controls['ref_date'] = controls['ref_date'].dt.normalize()   # [FIX-11]
controls['grade'] = 0

patient_df = pd.concat([cases[['mrn', 'toxicity', 'grade', 'ref_date']],
                        controls[['mrn', 'toxicity', 'grade', 'ref_date']]],
                       ignore_index=True)

n0 = len(patient_df)
patient_df = patient_df[patient_df['ref_date'].notna()].copy()
record_drop("missing reference date", n0, len(patient_df))    # [FIX-8]

# --- ICI status, identical rule for every grade [FIX-6] --------------------
# ref_date inside an ICI LOT's [lot_start, lot_end_eff]. With lot_end_pad this
# captures the 32% of AEs that present after the LOT formally closes.
ici_lots = skeleton[skeleton['ici_status']][['mrn', 'lot_start', 'lot_end_eff']]
chk = patient_df.merge(ici_lots, on='mrn', how='left')
chk['overlap'] = ((chk['lot_start'] <= chk['ref_date']) &
                  (chk['ref_date'] <= chk['lot_end_eff'])).fillna(False)
ici_at_ref = (chk.groupby(['mrn', 'toxicity'])['overlap'].any()
                 .rename('ici_at_ref').reset_index())
ici_ever = skeleton.groupby('mrn')['ici_status'].any().rename('ici_ever')

patient_df = (patient_df.merge(ici_at_ref, on=['mrn', 'toxicity'], how='left')
                        .merge(ici_ever, on='mrn', how='left'))
patient_df[['ici_at_ref', 'ici_ever']] = \
    patient_df[['ici_at_ref', 'ici_ever']].fillna(False)

print(f"\n  Patient-toxicity records: {len(patient_df):,}")
disagree = patient_df['ici_ever'] != patient_df['ici_at_ref']
print(f"  ici_ever vs ici_at_ref disagree: {disagree.sum():,} "
      f"({100 * disagree.mean():.1f}%)")
print(patient_df.groupby('grade')[['ici_ever', 'ici_at_ref']].mean().round(3))


# ============================================================================
# 3. LOAD MEDICATIONS
# ============================================================================

print("\n" + "=" * 74)
print("3. LOADING MEDICATIONS")
print("=" * 74)

med_frames = []

print("\n  rx_river.csv")
rx = pd.read_csv(DATA / "rx_river.csv", encoding='latin-1', low_memory=False)
n0 = len(rx)
rx['MRN'] = rx['MRN'].apply(standardize_mrn)
rx = rx[rx['MRN'].notna()].copy()
record_drop("rx: unparseable MRN", n0, len(rx))

rx['med_name'] = rx['GENERIC_NAME'].fillna(rx['DRUG_NAME'])
rx['steroid'] = rx['med_name'].apply(identify_steroid)
n0 = len(rx)
rx = rx[rx['steroid'].notna()].copy()
record_drop("rx: not a systemic steroid", n0, len(rx))

# [FIX-11] normalize
rx['start'] = pd.to_datetime(rx['START_DATE'], errors='coerce').dt.normalize()
rx['end_incl'] = pd.to_datetime(rx['END_DT'].replace('(null)', np.nan),
                                errors='coerce').dt.normalize()

n0 = len(rx)
rx = rx[rx['start'].notna()].copy()
record_drop("rx: missing start date", n0, len(rx))

# [FIX-12] start + 6 days is a 7-day inclusive course. v1 used 7 -> 8 days.
n_missing_end = int(rx['end_incl'].isna().sum())
rx['end_incl'] = rx['end_incl'].fillna(rx['start'] + timedelta(days=6))
print(f"    rx: imputed 7-day course for {n_missing_end:,} rows missing END_DT")

n_bad = int((rx['end_incl'] < rx['start']).sum())
if n_bad:
    print(f"    rx: {n_bad:,} rows with end < start — clamped to single day")
    rx.loc[rx['end_incl'] < rx['start'], 'end_incl'] = rx['start']

# [FIX-2] inclusive last day -> half-open upper bound
rx['end_ex'] = rx['end_incl'] + timedelta(days=1)

n_missing_dose = int(pd.to_numeric(rx['DOSE_AMT'], errors='coerce').isna().sum())
rx['dose'] = pd.to_numeric(rx['DOSE_AMT'], errors='coerce').fillna(DEFAULT_DOSE_MG)
print(f"    rx: imputed {DEFAULT_DOSE_MG:.0f}mg for {n_missing_dose:,} rows")

# [FIX-4] dosing frequency
freq_col = RX_FREQ_COL or next((c for c in RX_FREQ_CANDIDATES if c in rx.columns), None)
if freq_col:
    rx['doses_per_day'] = rx[freq_col].apply(parse_doses_per_day)
    print(f"    rx: frequency from '{freq_col}'; "
          f"{100 * (rx['doses_per_day'] == 1.0).mean():.0f}% defaulted to once-daily")
    print(f"    rx: {rx['doses_per_day'].value_counts().head(6).to_dict()}")
else:
    rx['doses_per_day'] = 1.0
    print("    rx: *** NO FREQUENCY COLUMN FOUND — assuming once-daily ***")
    print(f"    rx: columns available -> {sorted(rx.columns.tolist())}")

rx['mg_pred_per_day'] = (rx['dose'] * rx['doses_per_day']
                         * rx['steroid'].map(STEROID_EQUIVALENTS))
rx['source'] = 'rx'
med_frames.append(rx[['MRN', 'steroid', 'start', 'end_ex', 'mg_pred_per_day', 'source']])

print("\n  emar_river.csv")
emar = pd.read_csv(DATA / "emar_river.csv", encoding='latin-1', low_memory=False)
n0 = len(emar)
emar['MRN'] = emar['MRN'].apply(standardize_mrn)
emar = emar[emar['MRN'].notna()].copy()
record_drop("emar: unparseable MRN", n0, len(emar))

emar['steroid'] = emar['OO_ORD_NAME'].apply(identify_steroid)
n0 = len(emar)
emar = emar[emar['steroid'].notna()].copy()
record_drop("emar: not a systemic steroid", n0, len(emar))

emar['start'] = pd.to_datetime(emar['OOTO_SIGNIFIC_DTE'],
                               errors='coerce').dt.normalize()   # [FIX-11]
n0 = len(emar)
emar = emar[emar['start'].notna()].copy()
record_drop("emar: missing administration date", n0, len(emar))

emar['end_ex'] = emar['start'] + timedelta(days=1)   # [FIX-2] exactly one day

n_missing_dose = int(pd.to_numeric(emar['OOTO_TASK_DOSE'], errors='coerce').isna().sum())
emar['dose'] = pd.to_numeric(emar['OOTO_TASK_DOSE'],
                             errors='coerce').fillna(DEFAULT_DOSE_MG)
print(f"    emar: imputed {DEFAULT_DOSE_MG:.0f}mg for {n_missing_dose:,} rows")

# eMAR rows are individual administrations; frequency is implicit in the rows.
emar['mg_pred_per_day'] = emar['dose'] * emar['steroid'].map(STEROID_EQUIVALENTS)
emar['source'] = 'emar'
med_frames.append(emar[['MRN', 'steroid', 'start', 'end_ex', 'mg_pred_per_day', 'source']])

med_df = pd.concat(med_frames, ignore_index=True)
print(f"\n  Total steroid records: {len(med_df):,} "
      f"({(med_df['source'] == 'rx').sum():,} rx / "
      f"{(med_df['source'] == 'emar').sum():,} emar)")
print("\n  Records by steroid:")
for name, n in med_df['steroid'].value_counts().items():
    print(f"    {name:<20} {n:>9,}  (x{STEROID_EQUIVALENTS[name]:.2f} pred equiv)")


# ============================================================================
# 4. DAILY EXPOSURE TABLE, DEDUPLICATED ACROSS SOURCES
# ============================================================================

print("\n" + "=" * 74)
print("4. BUILDING DAILY EXPOSURE TABLE")
print("=" * 74)

bounds = (patient_df.groupby('mrn')['ref_date']
          .agg(win_lo='min', win_hi='max'))
bounds['win_lo'] -= timedelta(days=WINDOW_DAYS)
bounds['win_hi'] += timedelta(days=WINDOW_DAYS)

n0 = len(med_df)
med_df = med_df.merge(bounds, left_on='MRN', right_index=True, how='inner')
record_drop("meds for patients not in cohort", n0, len(med_df))

n0 = len(med_df)
med_df = med_df[(med_df['end_ex'] > med_df['win_lo']) &
                (med_df['start'] < med_df['win_hi'])].copy()
record_drop("meds outside any analysis window", n0, len(med_df))

med_df['start'] = med_df[['start', 'win_lo']].max(axis=1)
med_df['end_ex'] = med_df[['end_ex', 'win_hi']].min(axis=1)

daily = expand_to_days(med_df)
print(f"  Patient-steroid-days before dedup: {len(daily):,}")

# Multiple eMAR administrations of the same drug on one day are real: sum them.
daily = (daily.groupby(['MRN', 'steroid', 'day', 'source'], as_index=False)['mg_pred']
              .sum())

# [FIX-3] Collapse across sources: an inpatient course appears as both a
# prescription and administrations. Keep the eMAR value.
daily['src_rank'] = (daily['source'] == 'emar').astype(int)
before = len(daily)
daily = (daily.sort_values(['MRN', 'steroid', 'day', 'src_rank'])
              .drop_duplicates(subset=['MRN', 'steroid', 'day'], keep='last'))
print(f"  After cross-source dedup: {len(daily):,} "
      f"(removed {before - len(daily):,} rx days covered by eMAR)")

# Different steroids on the same day do legitimately add.
daily_total = (daily.groupby(['MRN', 'day'], as_index=False)['mg_pred'].sum()
                    .rename(columns={'mg_pred': 'mg'}))
print(f"  Patient-days with any steroid: {len(daily_total):,}")


# ============================================================================
# 5. WINDOW SUMS AND RATIO
# ============================================================================

print("\n" + "=" * 74)
print("5. POST/PRE RATIO")
print("=" * 74)
print(f"  Windows: pre = [ref-{WINDOW_DAYS}, ref)   post = [ref, ref+{WINDOW_DAYS})")
print("  Half-open and disjoint: the reference day belongs to post only. [FIX-10]")

merged = patient_df.merge(daily_total, left_on='mrn', right_on='MRN', how='left')
delta = (merged['day'] - merged['ref_date']).dt.days

merged['pre_mg'] = np.where(delta.between(-WINDOW_DAYS, -1), merged['mg'], 0.0)
merged['post_mg'] = np.where(delta.between(0, WINDOW_DAYS - 1), merged['mg'], 0.0)

agg = (merged.groupby(['mrn', 'toxicity'], as_index=False)[['pre_mg', 'post_mg']].sum())

results_df = patient_df.merge(agg, on=['mrn', 'toxicity'], how='left')
results_df[['pre_mg', 'post_mg']] = results_df[['pre_mg', 'post_mg']].fillna(0.0)
results_df = results_df.rename(columns={'pre_mg': 'pre_steroid',
                                        'post_mg': 'post_steroid'})

# [FIX-7] Flag no-exposure patients; do not park them at 0.0.
results_df['pre_any'] = results_df['pre_steroid'] > 0
results_df['post_any'] = results_df['post_steroid'] > 0
results_df['any_steroid'] = results_df['pre_any'] | results_df['post_any']
results_df['log2_ratio'] = np.where(
    results_df['any_steroid'],
    np.log2((results_df['post_steroid'] + EPSILON) /
            (results_df['pre_steroid'] + EPSILON)),
    np.nan)

print(f"\n  Total results: {len(results_df):,}")
print(f"  With any steroid: {results_df['any_steroid'].sum():,} "
      f"({100 * results_df['any_steroid'].mean():.1f}%)")

# Do not pickle a pandas DataFrame: pandas 3 StringDtype / DatetimeArray
# cannot be loaded by older Jupyter kernels. Store numpy arrays only.
def dump_cache(df, path):
    arrays, kinds = {}, {}
    for col in df.columns:
        s = df[col]
        if pd.api.types.is_datetime64_any_dtype(s):
            arrays[col] = np.asarray(s.astype("datetime64[ns]")).view("int64")
            kinds[col] = "datetime_ns"
        elif pd.api.types.is_bool_dtype(s):
            arrays[col] = np.asarray(s, dtype=bool)
            kinds[col] = "bool"
        elif pd.api.types.is_integer_dtype(s):
            arrays[col] = np.asarray(s, dtype=np.int64)
            kinds[col] = "int"
        elif pd.api.types.is_float_dtype(s):
            arrays[col] = np.asarray(s, dtype=np.float64)
            kinds[col] = "float"
        else:
            arrays[col] = np.asarray(s.astype(str), dtype=object)
            kinds[col] = "object"
    payload = {
        "_compat": "numpy_columns_v1",
        "columns": list(df.columns),
        "arrays": arrays,
        "kinds": kinds,
    }
    with open(path, "wb") as f:
        pickle.dump(payload, f, protocol=4)


out_path = DATA / "steroid_ratio_cache_LOT.pkl"
dump_cache(results_df, out_path)
print(f"  Saved: {out_path.name}")


# ============================================================================
# 6. VERIFICATION
# ============================================================================

print("\n" + "=" * 74)
print("6. COHORT COUNTS — exclusive (==k) and nested (>=k)")
print("=" * 74)
for tox in TOXICITIES:
    sub = results_df[results_df['toxicity'] == tox]
    print(f"\n{tox.upper()}:")
    for g in sorted(sub['grade'].unique()):
        excl = sub[sub['grade'] == g]
        nest = sub[sub['grade'] >= g]
        print(f"  Grade =={g}: n={len(excl):>7,} "
              f"({100 * excl['ici_at_ref'].mean():>3.0f}% ICI+ at ref) | "
              f">={g}: n={len(nest):>7,}")

print("\n" + "=" * 74)
print("7. ICI+ WITH STEROID EXPOSURE")
print("=" * 74)
ici_df = results_df[results_df['ici_at_ref'] & results_df['any_steroid']]
for tox in TOXICITIES:
    sub = ici_df[ici_df['toxicity'] == tox]
    print(f"\n{tox.upper()}: {len(sub):,} ICI+ patients with steroid")
    for g in sorted(sub['grade'].unique()):
        gd = sub[sub['grade'] == g]
        print(f"  Grade {g}: n={len(gd):>6,}  "
              f"pre={gd['pre_steroid'].median():>8.1f}  "
              f"post={gd['post_steroid'].median():>8.1f} mg-days  "
              f"log2={gd['log2_ratio'].median():>6.2f}")

print("\n" + "=" * 74)
print("8. SENSITIVITY — old ICI definition (ici_ever)")
print("=" * 74)
old_df = results_df[results_df['ici_ever'] & results_df['any_steroid']]
for tox in TOXICITIES:
    print(f"  {tox}: ici_ever n={len(old_df[old_df['toxicity'] == tox]):,} "
          f"vs ici_at_ref n={len(ici_df[ici_df['toxicity'] == tox]):,}")

print("\n" + "=" * 74)
print("DROP SUMMARY")
print("=" * 74)
for step, n in DROPS.items():
    print(f"  {step}: {n:,}")