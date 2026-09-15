#!/usr/bin/env python3
"""
Generate lab_results_cache.pkl for Figure 1D, E, F violin plots.
Run this script from terminal, then use the notebook for plotting.

Usage:
    python generate_lab_cache.py
"""
import os
import pickle
from pathlib import Path
from datetime import timedelta
import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

# ---------------------------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).parent.resolve()
ROOT = SCRIPT_DIR.parent
FIGURES = ROOT.parent.parent
DATA = FIGURES / "figures_data" / "figure 1" / "data"

GRADE_PATH = DATA / "grade_results_84k_FIXED_FP.csv"
LAB_PATH = DATA / "lab_results.csv"
CACHE_PATH = DATA / "lab_results_cache.pkl"

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------
TOXICITY_COLUMNS = ['pneumonitis', 'adrenal_insufficiency', 'liver_toxicity', 
                    'colitis', 'hyperthyroidism', 'hypothyroidism']

FILTERED_LABS = {
    'Cortisol': ['adrenal_insufficiency'],
    'ALT': ['liver_toxicity'],
    'AST': ['liver_toxicity'],
    'Alkaline_Phosphatase': ['liver_toxicity'],
    'Bilirubin_Total': ['liver_toxicity'],
    'TSH': ['hypothyroidism', 'hyperthyroidism'],
}

LAB_NAME_STANDARDIZATION = {
    'Alanine Aminotransferase (ALT), Plasma': 'ALT', 'Alanine Aminotransferase (ALT)': 'ALT',
    'Alanine Aminotransferase': 'ALT', 'ALT, Plasma': 'ALT', 'ALT': 'ALT', 'Alanine': 'ALT',
    'Aspartate Aminotransferase (AST), Plasma': 'AST', 'Aspartate Aminotransferase (AST)': 'AST',
    'Aspartate Aminotransferase': 'AST', 'AST, Plasma': 'AST', 'AST': 'AST',
    'Alkaline Phosphatase (ALK), Plasma': 'Alkaline_Phosphatase', 'Alkaline Phosphatase (ALK)': 'Alkaline_Phosphatase',
    'Alkaline Phosphatase, Plasma': 'Alkaline_Phosphatase', 'Alkaline Phosphatase': 'Alkaline_Phosphatase',
    'Bone Alkaline Phosphatase, S': None, 'Bone Alkaline Phosphatase': None,
    'Bilirubin, Total Plasma': 'Bilirubin_Total', 'Bilirubin, Total': 'Bilirubin_Total',
    'Total Bilirubin, Quest': 'Bilirubin_Total',
    'Thyroid Stimulating Hormone (TSH)': 'TSH', 'TSH': 'TSH',
    'Cortisol Level': 'Cortisol', 'Cortisol, Saliva': 'Cortisol', 'Cortisol, Free, 24 Hour, Urine': 'Cortisol',
}

GRADE_COLUMN_MAPPING = {
    'adrenal insufficiency': 'adrenal_insufficiency',
    'liver toxicity': 'liver_toxicity',
}

# ---------------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------------
def standardize_mrn(mrn):
    if pd.isna(mrn):
        return None
    mrn_str = ''.join(filter(str.isdigit, str(mrn).strip().replace('"', '').replace("'", '')))
    return mrn_str.zfill(8) if len(mrn_str) > 0 else None


def load_grade_predictions(grade_path):
    print("\nLoading grade predictions...")
    df = pd.read_csv(grade_path, encoding='latin-1')
    print(f"Loaded {len(df):,} records from {df['mrn'].nunique():,} patients")
    
    df['mrn'] = df['mrn'].apply(standardize_mrn)
    df = df[df['mrn'].notna()]
    df['window_start'] = pd.to_datetime(df['window_start'], errors='coerce')
    df['window_end'] = pd.to_datetime(df['window_end'], errors='coerce')
    df = df[df['window_start'].notna() & df['window_end'].notna()]
    
    for old_name, new_name in GRADE_COLUMN_MAPPING.items():
        if old_name in df.columns:
            df[new_name] = df[old_name]
    
    for tox in TOXICITY_COLUMNS:
        if tox in df.columns:
            df[tox] = pd.to_numeric(df[tox], errors='coerce').fillna(0).astype(int)
    
    print(f"After processing: {len(df):,} records from {df['mrn'].nunique():,} patients")
    return df


def classify_patients_fixed(llm_df_all):
    print(f"\n{'='*80}")
    print(f"CLASSIFYING PATIENTS - MAX GRADE ACROSS ALL LINES (FIXED)")
    print(f"{'='*80}")
    
    llm_df_all = llm_df_all.sort_values(['mrn', 'window_start'])
    
    results = []
    all_patient_mrns = set(llm_df_all['mrn'].unique())
    llm_by_patient = {mrn: group for mrn, group in llm_df_all.groupby('mrn')}
    
    for mrn in tqdm(all_patient_mrns, desc="Classifying patients"):
        patient_data = llm_by_patient[mrn]
        
        for tox in TOXICITY_COLUMNS:
            tox_grades = patient_data[tox].values
            max_grade = int(tox_grades.max())
            
            if max_grade > 0:
                first_max_idx = patient_data[patient_data[tox] == max_grade].index[0]
                first_max_window = patient_data.loc[first_max_idx]
                
                ae_date = first_max_window['window_start']
                post_start = ae_date
                post_end = ae_date + timedelta(days=60)
                pre_end = ae_date
                pre_start = ae_date - timedelta(days=60)
                
                results.append({
                    'mrn': mrn, 'ae_type': tox, 'grade': max_grade,
                    'ae_date': ae_date,
                    'post_start': post_start, 'post_end': post_end,
                    'pre_start': pre_start, 'pre_end': pre_end
                })
            else:
                all_starts = patient_data['window_start'].min()
                all_ends = patient_data['window_end'].max()
                midpoint = all_starts + (all_ends - all_starts) / 2
                
                post_start = midpoint
                post_end = midpoint + timedelta(days=60)
                pre_end = midpoint
                pre_start = midpoint - timedelta(days=60)
                
                results.append({
                    'mrn': mrn, 'ae_type': tox, 'grade': 0,
                    'ae_date': pd.NaT,
                    'post_start': post_start, 'post_end': post_end,
                    'pre_start': pre_start, 'pre_end': pre_end
                })
    
    df = pd.DataFrame(results)
    
    grade0_with_pre = df[(df['grade'] == 0) & df['pre_start'].notna()]
    print(f"\nFIX VERIFICATION: Grade 0 with pre_start defined: {len(grade0_with_pre):,}/{(df['grade'] == 0).sum():,}")
    print(f"Classified {len(df):,} patient-AE combinations from {len(all_patient_mrns):,} patients")
    
    return df


def collect_time_windows(patient_classifications):
    windows = []
    for _, row in patient_classifications.iterrows():
        if pd.notna(row['post_start']) and pd.notna(row['post_end']):
            windows.append({'mrn': row['mrn'], 'start': row['post_start'], 'end': row['post_end']})
        if pd.notna(row['pre_start']) and pd.notna(row['pre_end']):
            windows.append({'mrn': row['mrn'], 'start': row['pre_start'], 'end': row['pre_end']})
    print(f"Collected {len(windows):,} time windows")
    return windows


def load_lab_data(lab_path, patient_mrns, time_windows):
    print("\nLoading lab data (this may take a few minutes)...")
    
    mrn_windows = {}
    for window in time_windows:
        mrn = window['mrn']
        if mrn not in mrn_windows:
            mrn_windows[mrn] = []
        mrn_windows[mrn].append((window['start'], window['end']))
    
    df = pd.read_csv(lab_path, sep='\t', encoding='latin-1', low_memory=False, 
                     on_bad_lines='skip', quotechar='"')
    print(f"Initial load: {len(df):,} records")
    
    df['MRN'] = df['MRN'].apply(standardize_mrn)
    df = df[df['MRN'].notna() & df['MRN'].isin(patient_mrns)]
    print(f"After patient filter: {len(df):,} records")
    
    df['PERFORMED_DTE'] = pd.to_datetime(df['PERFORMED_DTE'], errors='coerce')
    df['RESULT_VALUE_NUMERIC'] = pd.to_numeric(df['RESULT_VALUE'], errors='coerce')
    df = df[df['PERFORMED_DTE'].notna() & df['RESULT_VALUE_NUMERIC'].notna()]
    
    print("Filtering to analysis time windows...")
    def is_in_window(row):
        mrn = row['MRN']
        date = row['PERFORMED_DTE']
        if mrn not in mrn_windows:
            return False
        for start, end in mrn_windows[mrn]:
            if start <= date <= end:
                return True
        return False
    
    df = df[df.apply(is_in_window, axis=1)]
    print(f"After time window filter: {len(df):,} records")
    
    df['STANDARDIZED_TEST'] = df['SUBTEST_NAME'].map(LAB_NAME_STANDARDIZATION)
    df = df[df['STANDARDIZED_TEST'].notna() & df['STANDARDIZED_TEST'].isin(FILTERED_LABS.keys())]
    print(f"After lab filter: {len(df):,} records ({df['STANDARDIZED_TEST'].nunique()} unique tests)")
    
    return df


def calculate_auc(mrn, lab_test, start_date, end_date, lab_df):
    if pd.isna(start_date) or pd.isna(end_date) or start_date >= end_date:
        return np.nan
    
    patient_labs = lab_df[
        (lab_df['MRN'] == mrn) &
        (lab_df['STANDARDIZED_TEST'] == lab_test) &
        (lab_df['PERFORMED_DTE'] >= start_date) &
        (lab_df['PERFORMED_DTE'] <= end_date)
    ].sort_values('PERFORMED_DTE')
    
    if len(patient_labs) < 1:
        return np.nan
    if len(patient_labs) == 1:
        return patient_labs['RESULT_VALUE_NUMERIC'].iloc[0]
    
    times = (patient_labs['PERFORMED_DTE'] - start_date).dt.total_seconds() / (24 * 3600)
    values = patient_labs['RESULT_VALUE_NUMERIC'].values
    auc = np.trapezoid(values, times)
    
    total_days = (end_date - start_date).total_seconds() / (24 * 3600)
    return auc / total_days if total_days > 0 else auc


def process_lab_test(args):
    lab_test, ae_patients, lab_subset, task_id = args
    results = []
    
    total = len(ae_patients)
    for i, (_, row) in enumerate(ae_patients.iterrows()):
        if i % 1000 == 0:
            print(f"  [{lab_test[:8]:8s}] {i:,}/{total:,} ({100*i/total:.0f}%)", flush=True)
        
        post_auc = calculate_auc(row['mrn'], lab_test, row['post_start'], row['post_end'], lab_subset)
        pre_auc = calculate_auc(row['mrn'], lab_test, row['pre_start'], row['pre_end'], lab_subset)
        
        # Keep anyone with a post-window AUC. Exact zeros are not expected for
        # these labs; log2_post uses +0.01 so a zero would still be defined.
        # (The old `post_auc > 0` cut was leftover from the ratio.)
        if pd.notna(post_auc):
            results.append({
                'mrn': row['mrn'], 'ae_type': row['ae_type'], 'lab_test': lab_test,
                'grade': row['grade'],
                'post_auc': post_auc,
                'pre_auc': pre_auc if pd.notna(pre_auc) and pre_auc > 0 else np.nan
            })
    
    print(f"  [{lab_test[:8]:8s}] DONE - {len(results):,} results", flush=True)
    return results


def analyze_lab_values(patient_classifications, lab_df, use_parallel=True):
    print(f"\n{'='*60}")
    print(f"CALCULATING LAB AUCs (POST AND PRE)")
    print(f"{'='*60}")
    
    tasks = []
    task_id = 0
    for lab_test in lab_df['STANDARDIZED_TEST'].unique():
        if lab_test not in FILTERED_LABS:
            continue
        lab_subset = lab_df[lab_df['STANDARDIZED_TEST'] == lab_test].copy()
        if len(lab_subset) == 0:
            continue
        
        for ae_type in FILTERED_LABS.get(lab_test, []):
            ae_patients = patient_classifications[patient_classifications['ae_type'] == ae_type]
            if len(ae_patients) > 0:
                print(f"  {lab_test} / {ae_type}: {len(ae_patients):,} patients")
                tasks.append((lab_test, ae_patients, lab_subset, task_id))
                task_id += 1
    
    print(f"\nProcessing {len(tasks)} lab test-AE combinations...")
    
    all_results = []
    if use_parallel and len(tasks) > 1:
        n_cores = min(cpu_count() - 1, len(tasks), 8)
        print(f"Using {n_cores} cores (parallel)\n")
        with Pool(processes=n_cores) as pool:
            for results_list in pool.imap_unordered(process_lab_test, tasks):
                all_results.extend(results_list)
    else:
        print("Processing serially...")
        for task in tasks:
            all_results.extend(process_lab_test(task))
    
    df = pd.DataFrame(all_results)
    print(f"\nCalculated AUCs: {len(df):,} records")
    
    df['log2_post'] = np.log2(df['post_auc'] + 0.01)
    df['log2_ratio'] = np.nan
    valid_ratio = df['pre_auc'].notna() & (df['pre_auc'] > 0)
    df.loc[valid_ratio, 'log2_ratio'] = np.log2(df.loc[valid_ratio, 'post_auc'] / df.loc[valid_ratio, 'pre_auc'])
    
    print(f"  Valid log2_post: {df['log2_post'].notna().sum():,}")
    print(f"  Valid log2_ratio: {df['log2_ratio'].notna().sum():,}")
    
    return df


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("="*80)
    print("GENERATING LAB RESULTS CACHE")
    print("="*80)
    
    print(f"\nPaths:")
    print(f"  Grade file: {GRADE_PATH}")
    print(f"  Lab file: {LAB_PATH}")
    print(f"  Output cache: {CACHE_PATH}")
    
    # Step 1: Load grade predictions
    llm_df_all = load_grade_predictions(GRADE_PATH)
    
    # Step 2: Classify patients
    patient_classifications = classify_patients_fixed(llm_df_all)
    
    # Step 3: Collect time windows
    time_windows = collect_time_windows(patient_classifications)
    
    # Step 4: Load and filter lab data
    patient_mrns = set(patient_classifications['mrn'].unique())
    lab_df = load_lab_data(LAB_PATH, patient_mrns, time_windows)
    
    # Step 5: Calculate AUCs (with multiprocessing)
    results_df = analyze_lab_values(patient_classifications, lab_df, use_parallel=True)
    
    # Step 6: Save cache
    with open(CACHE_PATH, 'wb') as f:
        pickle.dump(results_df, f)
    
    print(f"\n{'='*80}")
    print(f"DONE! Saved cache: {CACHE_PATH}")
    print(f"{'='*80}")
    
    # Verification
    print(f"\nVerification:")
    grade0 = results_df[results_df['grade'] == 0]
    print(f"  Grade 0 with valid log2_ratio: {grade0['log2_ratio'].notna().sum():,}/{len(grade0):,}")
    
    print(f"\n  By grade:")
    for grade in sorted(results_df['grade'].unique()):
        g = results_df[results_df['grade'] == grade]
        valid = g['log2_ratio'].notna().sum()
        print(f"    Grade {int(grade)}: {valid:,}/{len(g):,} ({100*valid/len(g):.1f}%)")


if __name__ == "__main__":
    main()
