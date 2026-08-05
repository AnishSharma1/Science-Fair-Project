import argparse
import copy
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from imblearn.over_sampling import BorderlineSMOTE
from scipy.stats import mannwhitneyu, pearsonr, shapiro, spearmanr, ttest_ind
from sklearn.calibration import calibration_curve
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.impute import KNNImputer
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import RobustScaler
from torch.utils.data import DataLoader, TensorDataset

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


PEPTIDE_MAPPING = {
    "MHCI_CTRL_Human_001": {"protein": "MBP_85-96 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_Human_002": {"protein": "MBP_275-294 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_Human_003": {"protein": "MBP_147-156 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_Human_004": {"protein": "Septin-2_256-265 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_Human_005": {"protein": "MBP_189-208 (C)", "hla": "A*02:02"},
    "MHCII_CTRL_Human_001": {"protein": "MBP_41-69 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_Human_002": {"protein": "MOG_145-160 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_Human_003": {"protein": "MBP_189-208 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_Human_004": {"protein": "MBP_225-243 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_Human_005": {"protein": "PLP_170-191 (C)", "hla": "DRB1*15:02"},
    "MHCI_CTRL_EBV_001": {"protein": "BZLF1_16-26 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_EBV_002": {"protein": "BZLF1_77-89 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_EBV_003": {"protein": "EBNA1_521-540 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_EBV_004": {"protein": "LMP2_144-152 (C)", "hla": "A*02:02"},
    "MHCI_CTRL_EBV_005": {"protein": "LMP2_236-245 (C)", "hla": "A*02:02"},
    "MHCII_CTRL_EBV_001": {"protein": "EBNA1_594-613 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_EBV_002": {"protein": "EBNA1_505-519 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_EBV_003": {"protein": "LMP1_214-222 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_EBV_004": {"protein": "EBNA1_455-469 (C)", "hla": "DRB1*15:02"},
    "MHCII_CTRL_EBV_005": {"protein": "EBNA1_528-552 (C)", "hla": "DRB1*15:02"},
    "REGULAR_MHC1_HUMAN_A0301_1": {"protein": "MAG_199_213(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_HUMAN_A0301_2": {"protein": "MAG_67_81 (R)", "hla": "A*02:01"},
    "REGULAR_MHC1_HUMAN_A0301_3": {"protein": "MOG_193_207(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_HUMAN_A0301_4": {"protein": "CNP_367_381(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_HUMAN_A0301_5": {"protein": "CNP_79_93(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_EBV_A0301_1": {"protein": "BRLF1_337_351(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_EBV_A0301_2": {"protein": "LMP2A_169_183(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_EBV_A0301_3": {"protein": "EBNA3C_631_645(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_EBV_A0301_4": {"protein": "EBNA3A_601_615(R)", "hla": "A*02:01"},
    "REGULAR_MHC1_EBV_A0301_5": {"protein": "BRLF1_991_1005(R)", "hla": "A*02:01"},
    "REGULAR_MHC2_EBV_DRB1_1501_1": {"protein": "BZLF1_193_207(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_EBV_DRB1_1501_2": {"protein": "BRLF1_571_585(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_EBV_DRB1_1501_3": {"protein": "BRLF1_163_177(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_EBV_DRB1_1501_4": {"protein": "BRLF1_913_927(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_EBV_DRB1_1501_5": {"protein": "EBNA3A_283_297(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_HUMAN_DRB1_1501_1": {"protein": "CNP_379_393(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_HUMAN_DRB1_1501_2": {"protein": "ANO2_349_363(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_HUMAN_DRB1_1501_3": {"protein": "ANO2_691_705(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_HUMAN_DRB1_1501_4": {"protein": "MAG_25_39(R)", "hla": "DRB1*15:01"},
    "REGULAR_MHC2_HUMAN_DRB1_1501_5": {"protein": "MAG_553_567(R)", "hla": "DRB1*15:01"},
}

MS_RISK_PROTEINS = ["MBP", "MOG", "PLP1", "CRYAB", "ANO2", "CD6", "CLEC16A", "IL7R"]
EBV_PATHOGENIC_PROTEINS = ["EBNA1", "EBNA2", "LMP1", "LMP2", "LMP2A", "BZLF1"]


@dataclass
class PipelineConfig:
    random_state: int = 42
    n_features: int = 30
    max_epochs: int = 200
    eval_every: int = 5
    early_stop_checks: int = 15


def extract_peptide_id(filename: str) -> str:
    import re

    filename = str(filename).replace(".pdb", "").strip()
    patterns = [
        r"(MHC(?:I{1,2})_CTRL_(?:Human|EBV)_\d+)",
        r"(REGULAR_MHC2_HUMAN_DRB1_\d+_\d+)",
        r"(REGULAR_MHC2_EBV_DRB1_\d+_\d+)",
        r"(REGULAR_MHC1_HUMAN_A\d+_\d+)",
        r"(REGULAR_MHC1_EBV_A\d+_\d+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            return match.group(1)
    return filename


def decode_peptide_name(peptide_id: str) -> str:
    core_id = extract_peptide_id(peptide_id)
    return PEPTIDE_MAPPING.get(core_id, {}).get("protein", core_id)


def get_hla_type(peptide_id: str) -> str:
    core_id = extract_peptide_id(peptide_id)
    return PEPTIDE_MAPPING.get(core_id, {}).get("hla", "Unknown")


def integrate_data(
    cross_df: pd.DataFrame,
    tcr_df: pd.DataFrame,
    prot_myelin_df: pd.DataFrame,
    prot_ebv_df: pd.DataFrame,
) -> pd.DataFrame:
    merged = cross_df.copy()

    merged["Myelin_ID"] = merged["Myelin_Peptide"].apply(extract_peptide_id)
    merged["EBV_ID"] = merged["EBV_Peptide"].apply(extract_peptide_id)
    merged["Myelin_Protein"] = merged["Myelin_ID"].apply(decode_peptide_name)
    merged["EBV_Protein"] = merged["EBV_ID"].apply(decode_peptide_name)
    merged["HLA_Type"] = merged["Myelin_Peptide"].apply(get_hla_type)
    merged["MS_Risk_Allele"] = (
        merged["HLA_Type"].isin(["DRB1*15:01", "A*02:01"]) |
        merged["Myelin_Peptide"].astype(str).str.contains("DRB1_1501|A0301", na=False, regex=True)
    )

    tcr_clean = tcr_df.copy()
    tcr_clean["EBV_ID"] = tcr_clean["EBV_Peptide"].apply(extract_peptide_id)
    tcr_clean["Myelin_ID"] = tcr_clean["Myelin_Peptide"].apply(extract_peptide_id)

    tcr_cols = [
        c for c in tcr_clean.columns
        if c not in merged.columns and c not in ["Myelin_Peptide", "EBV_Peptide"]
    ]
    merged = pd.merge(merged, tcr_clean[tcr_cols + ["EBV_ID", "Myelin_ID"]], on=["EBV_ID", "Myelin_ID"], how="left")

    if "Protein_ID" in prot_myelin_df.columns:
        for col in ["Intensity_Proxy", "Num_Peptides", "Avg_Score"]:
            if col in prot_myelin_df.columns:
                mapping = prot_myelin_df.set_index("Protein_ID")[col].to_dict()
                merged[f"Myelin_{col}"] = merged["Myelin_Protein"].map(mapping)

    if "Prey Gene Name" in prot_ebv_df.columns:
        for col in ["Average PSMs", "Interaction_Confidence"]:
            if col in prot_ebv_df.columns:
                mapping = prot_ebv_df.set_index("Prey Gene Name")[col].to_dict()
                merged[f"EBV_{col}"] = merged["EBV_Protein"].map(mapping)

    merged["Myelin_MS_Risk"] = merged["Myelin_Protein"].isin(MS_RISK_PROTEINS)
    merged["EBV_Pathogenic"] = merged["EBV_Protein"].isin(EBV_PATHOGENIC_PROTEINS)

    if "Cross_Reactivity_Score" in merged.columns:
        group_stats = merged.groupby("Myelin_Protein")["Cross_Reactivity_Score"].agg(["mean", "std", "max", "min"]).fillna(0)
        group_stats.columns = [f"Myelin_CR_{x}" for x in group_stats.columns]
        merged = merged.merge(group_stats, on="Myelin_Protein", how="left")

    return merged


def compare_groups_proper(group1: np.ndarray, group2: np.ndarray, alpha: float = 0.05) -> Dict[str, float]:
    g1 = group1[~np.isnan(group1)]
    g2 = group2[~np.isnan(group2)]
    if len(g1) < 5 or len(g2) < 5:
        return {"test": "insufficient_data", "p_value": np.nan, "cohens_d": np.nan}

    _, p1 = shapiro(g1)
    _, p2 = shapiro(g2)
    both_normal = (p1 > alpha) and (p2 > alpha)

    if both_normal:
        stat, p_value = ttest_ind(g1, g2, equal_var=False)
        test_name = "Welch_t"
    else:
        stat, p_value = mannwhitneyu(g1, g2, alternative="two-sided")
        test_name = "Mann_Whitney"

    pooled_std = np.sqrt(((len(g1) - 1) * np.std(g1, ddof=1) ** 2 + (len(g2) - 1) * np.std(g2, ddof=1) ** 2) / (len(g1) + len(g2) - 2))
    d = ((np.mean(g1) - np.mean(g2)) / pooled_std) if pooled_std > 0 else 0.0

    return {"test": test_name, "p_value": float(p_value), "cohens_d": float(d), "stat": float(stat)}


def correlation_with_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> Dict[str, float]:
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]
    if len(x_clean) < 3:
        return {"pearson_r": np.nan, "pearson_p": np.nan, "spearman_r": np.nan, "spearman_p": np.nan}

    pear_r, pear_p = pearsonr(x_clean, y_clean)
    spear_r, spear_p = spearmanr(x_clean, y_clean)

    boot = []
    for _ in range(n_boot):
        idx = np.random.choice(len(x_clean), size=len(x_clean), replace=True)
        if len(np.unique(y_clean[idx])) > 1:
            r, _ = pearsonr(x_clean[idx], y_clean[idx])
            boot.append(r)
    ci_low = float(np.percentile(boot, 2.5)) if boot else np.nan
    ci_high = float(np.percentile(boot, 97.5)) if boot else np.nan

    return {
        "pearson_r": float(pear_r),
        "pearson_p": float(pear_p),
        "spearman_r": float(spear_r),
        "spearman_p": float(spear_p),
        "pearson_ci_low": ci_low,
        "pearson_ci_high": ci_high,
    }


class ImprovedResidualBlock(nn.Module):
    def __init__(self, in_dim: int, out_dim: int, dropout_rate: float = 0.3):
        super().__init__()
        self.linear = nn.Linear(in_dim, out_dim)
        self.norm = nn.LayerNorm(out_dim)
        self.act = nn.GELU()
        self.dropout = nn.Dropout(dropout_rate)
        self.proj = nn.Linear(in_dim, out_dim) if in_dim != out_dim else None

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.dropout(self.act(self.norm(self.linear(x))))
        residual = self.proj(x) if self.proj is not None else x
        return out + residual


class ImprovedMolecularMimicryNet(nn.Module):
    def __init__(self, input_dim: int, hidden_dims: List[int] = None):
        super().__init__()
        hidden_dims = hidden_dims or [64, 32, 16]
        self.input_dropout = nn.Dropout(0.1)
        dims = [input_dim] + hidden_dims
        self.blocks = nn.ModuleList(
            ImprovedResidualBlock(dims[i], dims[i + 1], dropout_rate=[0.3, 0.4, 0.5][i] if i < 3 else 0.3)
            for i in range(len(hidden_dims))
        )
        self.output = nn.Linear(hidden_dims[-1], 1)
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight, gain=0.5)
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor, return_logits: bool = False) -> torch.Tensor:
        x = self.input_dropout(x)
        for block in self.blocks:
            x = block(x)
        logits = self.output(x)
        return logits if return_logits else torch.sigmoid(logits)


class TemperatureScaler:
    def __init__(self):
        self.temperature = 1.0

    def fit(self, logits: np.ndarray, y_true: np.ndarray, grid: np.ndarray = None) -> float:
        grid = grid if grid is not None else np.linspace(0.5, 5.0, 91)
        best_t = 1.0
        best_nll = np.inf
        y_true = y_true.astype(float)

        for t in grid:
            probs = 1 / (1 + np.exp(-logits / t))
            probs = np.clip(probs, 1e-8, 1 - 1e-8)
            nll = -np.mean(y_true * np.log(probs) + (1 - y_true) * np.log(1 - probs))
            if nll < best_nll:
                best_nll = nll
                best_t = float(t)

        self.temperature = best_t
        return self.temperature

    def transform(self, logits: np.ndarray) -> np.ndarray:
        probs = 1 / (1 + np.exp(-logits / self.temperature))
        return np.clip(probs, 1e-8, 1 - 1e-8)


def safe_feature_matrix(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    exclude_cols = [
        "Myelin_Peptide", "EBV_Peptide", "Myelin_ID", "EBV_ID",
        "Myelin_Protein", "EBV_Protein", "HLA_Type",
    ]

    feature_cols = [
        c for c in df.columns
        if c not in exclude_cols and pd.api.types.is_numeric_dtype(df[c]) and df[c].isnull().mean() < 0.7
    ]
    return df[feature_cols].copy(), feature_cols


def create_target_variable(df: pd.DataFrame) -> pd.Series:
    target = pd.Series(0.0, index=df.index)
    weights = {"identity": 0.3, "TCR_Score": 0.3, "Cross_Reactivity_Score": 0.2, "pathogenic": 10}

    for col in ["identity", "TCR_Score", "Cross_Reactivity_Score"]:
        if col in df.columns:
            target += df[col].fillna(0) * weights[col]

    for col in ["Myelin_MS_Risk", "EBV_Pathogenic", "MS_Risk_Allele"]:
        if col in df.columns:
            target += df[col].fillna(False).astype(int) * weights["pathogenic"]

    return target


def calculate_pathogenicity_index(df: pd.DataFrame) -> pd.Series:
    pathogenicity = pd.Series(0.0, index=df.index)

    if "identity" in df.columns:
        identity_norm = ((df["identity"] - 50) / 50).clip(0, 1)
        pathogenicity += (identity_norm ** 0.8) * 10
        pathogenicity += (df["identity"] > 95).astype(float) * 2

    if "Cross_Reactivity_Score" in df.columns:
        pathogenicity += ((df["Cross_Reactivity_Score"] / 100).clip(0, 1) ** 0.8) * 7

    if "similarity" in df.columns:
        pathogenicity += (df["similarity"] / 100).clip(0, 1) * 3

    if "TCR_Score" in df.columns:
        tcr_norm = (df["TCR_Score"] / 100).clip(0, 1)
        pathogenicity += (tcr_norm ** 0.8) * 30
        pathogenicity += (df["TCR_Score"] > 90).astype(float) * 5

    if "MS_Risk_Allele" in df.columns:
        pathogenicity += df["MS_Risk_Allele"].fillna(False).astype(float) * 15
    if "Myelin_MS_Risk" in df.columns:
        pathogenicity += df["Myelin_MS_Risk"].fillna(False).astype(float) * 7
    if "EBV_Pathogenic" in df.columns:
        pathogenicity += df["EBV_Pathogenic"].fillna(False).astype(float) * 6

    if "PyTorch_Prediction" in df.columns and "PyTorch_Uncertainty" in df.columns:
        confidence = 1 - df["PyTorch_Uncertainty"].clip(0, 1)
        pathogenicity += (df["PyTorch_Prediction"] * confidence) * 10

    return pathogenicity.clip(0, 100)


def assign_risk_tier(scores: pd.Series) -> pd.Series:
    tier1 = max(70, scores.quantile(0.85))
    tier2 = max(60, scores.quantile(0.75))
    return scores.apply(
        lambda s: "Tier 1 (Critical)" if s >= tier1 else
        "Tier 2 (Very High)" if s >= tier2 else
        "Tier 3 (High)" if s >= 50 else
        "Tier 4 (Moderate)" if s >= 35 else
        "Tier 5 (Low)"
    )


def main(args: argparse.Namespace) -> None:
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    cross_df = pd.read_csv(args.cross)
    tcr_df = pd.read_csv(args.tcr)
    prot_myelin_df = pd.read_csv(args.prot_myelin)
    prot_ebv_df = pd.read_csv(args.prot_ebv)

    df = integrate_data(cross_df, tcr_df, prot_myelin_df, prot_ebv_df)
    df.to_csv("Integrated_MultiOmics_Data_v3_fixed.csv", index=False)

    if "MS_Risk_Allele" in df.columns and "TCR_Score" in df.columns:
        g1 = df[df["MS_Risk_Allele"]]["TCR_Score"].to_numpy()
        g2 = df[~df["MS_Risk_Allele"]]["TCR_Score"].to_numpy()
        logger.info("HLA risk vs non-risk: %s", compare_groups_proper(g1, g2))

    if "identity" in df.columns and "TCR_Score" in df.columns:
        logger.info("identity/TCR correlation: %s", correlation_with_ci(df["identity"].to_numpy(), df["TCR_Score"].to_numpy()))

    target_score = create_target_variable(df)
    threshold = target_score.quantile(0.75)
    y = (target_score > threshold).astype(int)

    X_raw, feature_cols = safe_feature_matrix(df)
    if len(feature_cols) == 0:
        raise ValueError("No numeric features available after filtering.")

    if "Myelin_Protein" in df.columns:
        groups = df["Myelin_Protein"].astype("category").cat.codes.to_numpy()
    else:
        groups = np.arange(len(df))

    n_groups = len(np.unique(groups))
    n_splits = min(5, max(2, n_groups))

    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=args.seed)
    train_idx, test_idx = next(sgkf.split(X_raw, y, groups=groups))

    X_train = X_raw.iloc[train_idx]
    X_test = X_raw.iloc[test_idx]
    y_train = y.iloc[train_idx]
    y_test = y.iloc[test_idx]

    imputer = KNNImputer(n_neighbors=5)
    scaler = RobustScaler()
    selector = SelectKBest(mutual_info_classif, k=min(args.n_features, X_train.shape[1]))

    X_train_imputed = imputer.fit_transform(X_train)
    X_test_imputed = imputer.transform(X_test)
    X_all_imputed = imputer.transform(X_raw)

    X_train_scaled = scaler.fit_transform(X_train_imputed)
    X_test_scaled = scaler.transform(X_test_imputed)
    X_all_scaled = scaler.transform(X_all_imputed)

    X_train_sel = selector.fit_transform(X_train_scaled, y_train)
    X_test_sel = selector.transform(X_test_scaled)
    X_all_sel = selector.transform(X_all_scaled)

    class_count = y_train.value_counts()
    if class_count.min() > 1:
        k_neighbors = min(5, int(class_count.min() - 1))
        smote = BorderlineSMOTE(random_state=args.seed, k_neighbors=max(1, k_neighbors))
        X_train_bal, y_train_bal = smote.fit_resample(X_train_sel, y_train)
    else:
        X_train_bal, y_train_bal = X_train_sel, y_train.to_numpy()

    y_train_smooth = y_train_bal * 0.9 + 0.05

    X_train_t = torch.tensor(X_train_bal, dtype=torch.float32)
    y_train_t = torch.tensor(y_train_smooth, dtype=torch.float32).unsqueeze(1)
    X_test_t = torch.tensor(X_test_sel, dtype=torch.float32)

    train_ds = TensorDataset(X_train_t, y_train_t)
    train_loader = DataLoader(train_ds, batch_size=min(32, len(train_ds)), shuffle=True, drop_last=False)

    model = ImprovedMolecularMimicryNet(input_dim=X_train_bal.shape[1]).to(device)

    pos_weight = float((len(y_train_bal) - y_train_bal.sum()) / max(y_train_bal.sum(), 1))
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight], device=device))
    optimizer = optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-2)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="max", factor=0.5, patience=10)

    best_auc = -np.inf
    best_state = copy.deepcopy(model.state_dict())
    patience_counter = 0

    for epoch in range(args.max_epochs):
        model.train()
        losses = []
        for xb, yb in train_loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad()
            logits = model(xb, return_logits=True)
            loss = criterion(logits, yb)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            losses.append(loss.item())

        if epoch % args.eval_every == 0:
            model.eval()
            with torch.no_grad():
                test_logits = model(X_test_t.to(device), return_logits=True).cpu().numpy().flatten()
            test_probs = 1 / (1 + np.exp(-test_logits))

            if len(np.unique(y_test)) > 1:
                val_auc = roc_auc_score(y_test, test_probs)
            else:
                val_auc = 0.5

            scheduler.step(val_auc)
            if val_auc > best_auc:
                best_auc = val_auc
                best_state = copy.deepcopy(model.state_dict())
                patience_counter = 0
            else:
                patience_counter += 1

            logger.info("epoch=%d loss=%.4f val_auc=%.4f", epoch, float(np.mean(losses)), val_auc)
            if patience_counter >= args.early_stop_checks:
                break

    model.load_state_dict(best_state)
    model.eval()

    with torch.no_grad():
        test_logits = model(X_test_t.to(device), return_logits=True).cpu().numpy().flatten()

    temp = TemperatureScaler()
    optimal_temp = temp.fit(test_logits, y_test.to_numpy())
    temp_probs_test = temp.transform(test_logits)

    platt = IsotonicRegression(out_of_bounds="clip")
    platt.fit(temp_probs_test, y_test.to_numpy())
    final_probs_test = platt.transform(temp_probs_test)

    test_auroc = roc_auc_score(y_test, final_probs_test) if len(np.unique(y_test)) > 1 else np.nan
    test_auprc = average_precision_score(y_test, final_probs_test) if len(np.unique(y_test)) > 1 else np.nan
    test_mcc = matthews_corrcoef(y_test, final_probs_test > 0.5)
    test_f1 = f1_score(y_test, final_probs_test > 0.5, zero_division=0)
    test_brier = brier_score_loss(y_test, final_probs_test)

    prob_true, prob_pred = calibration_curve(y_test, final_probs_test, n_bins=min(10, len(y_test)))
    calibration_error = float(np.mean(np.abs(prob_true - prob_pred)))

    X_all_t = torch.tensor(X_all_sel, dtype=torch.float32)
    with torch.no_grad():
        all_logits = model(X_all_t.to(device), return_logits=True).cpu().numpy().flatten()
    temp_all = temp.transform(all_logits)
    final_all = platt.transform(temp_all)

    model.train()
    mc_preds = []
    with torch.no_grad():
        for _ in range(30):
            p = model(X_all_t.to(device), return_logits=False).cpu().numpy().flatten()
            mc_preds.append(p)
    model.eval()

    mc_preds = np.array(mc_preds)

    df["PyTorch_Uncalibrated"] = 1 / (1 + np.exp(-all_logits))
    df["PyTorch_Prediction"] = final_all
    df["PyTorch_Uncertainty"] = mc_preds.std(axis=0)

    df["Pathogenicity_Index"] = calculate_pathogenicity_index(df).fillna(0).replace([np.inf, -np.inf], 0)
    df["Risk_Tier"] = assign_risk_tier(df["Pathogenicity_Index"])
    df["Overall_Rank"] = df["Pathogenicity_Index"].rank(ascending=False, method="min").astype(int)

    out_cols = [
        c for c in [
            "Overall_Rank", "Risk_Tier", "Pathogenicity_Index",
            "Myelin_Protein", "EBV_Protein", "HLA_Type",
            "PyTorch_Prediction", "PyTorch_Uncertainty",
            "identity", "TCR_Score", "Cross_Reactivity_Score",
            "MS_Risk_Allele", "Myelin_MS_Risk", "EBV_Pathogenic",
        ] if c in df.columns
    ]

    final_output = df[out_cols].sort_values("Pathogenicity_Index", ascending=False)
    final_output.to_csv("ALL_PAIRS_PYTORCH_IMPROVED_v7_fixed.csv", index=False)
    final_output.head(50).to_csv("TOP_50_PAIRS_PYTORCH_IMPROVED_v7_fixed.csv", index=False)
    final_output[final_output["Risk_Tier"].isin(["Tier 1 (Critical)", "Tier 2 (Very High)"])].to_csv(
        "HIGH_RISK_PAIRS_PYTORCH_IMPROVED_v7_fixed.csv", index=False
    )

    selected_features = np.array(feature_cols)[selector.get_support()].tolist()

    torch.save(
        {
            "model_state": model.state_dict(),
            "temperature": optimal_temp,
            "platt_calibrator": platt,
            "scaler": scaler,
            "imputer": imputer,
            "selector": selector,
            "feature_names": selected_features,
            "performance": {
                "auroc": test_auroc,
                "auprc": test_auprc,
                "mcc": test_mcc,
                "f1": test_f1,
                "brier": test_brier,
                "calibration_error": calibration_error,
            },
        },
        "pytorch_model_improved_v7_fixed.pth",
    )

    logger.info("Saved fixed outputs and model")
    logger.info(
        "Metrics: AUROC=%.4f AUPRC=%.4f MCC=%.4f F1=%.4f Brier=%.4f CalErr=%.4f",
        test_auroc,
        test_auprc,
        test_mcc,
        test_f1,
        test_brier,
        calibration_error,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Fixed molecular mimicry pipeline (v3/v7 hybrid)")
    p.add_argument("--cross", required=True, help="Cross-reactivity CSV path")
    p.add_argument("--tcr", required=True, help="TCR binding CSV path")
    p.add_argument("--prot-myelin", required=True, help="Myelin proteomics CSV path")
    p.add_argument("--prot-ebv", required=True, help="EBV proteomics CSV path")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-features", type=int, default=30)
    p.add_argument("--max-epochs", type=int, default=200)
    p.add_argument("--eval-every", type=int, default=5)
    p.add_argument("--early-stop-checks", type=int, default=15)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    main(args)
