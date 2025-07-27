import os, re, logging
import pandas as pd
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions as Reactions

# --- Logging / RDKit noise ---
RDLogger.DisableLog('rdApp.*')
logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
log = logging.getLogger(__name__)

# --- Config ---
COLS = ['ReactionSMILES','PatentNumber','ParagraphNum','Year','TextMinedYield','CalculatedYield']
BAD_TOKENS_RE = re.compile(r'(\[\*|\[R|;|\?|\$\()')  # common SMARTS/dummy cues that explode parsers
CHUNKSIZE = 100_000

def safe_rxn_from_string(rsmi: str):
    """Parse USPTO reaction strings safely. Returns RDKit ChemicalReaction or None."""
    if not isinstance(rsmi, str):
        return None
    s = rsmi.strip()
    if not s or '>' not in s or BAD_TOKENS_RE.search(s):
        return None
    # 1) Correct parser for reaction SMILES (handles 'reactants>agents>products' and 'reactants>>products')
    try:
        rxn = Reactions.ReactionFromSmiles(s, useSmiles=True)
        if rxn and rxn.GetNumReactantTemplates() > 0 and rxn.GetNumProductTemplates() > 0:
            return rxn
    except Exception:
        pass
    # 2) Fallback: build from components
    try:
        parts = s.split('>')
        if len(parts) == 3:   # reactants>agents>products
            r_str, _, p_str = parts
        else:                  # reactants>>products
            r_str, p_str = s.split('>>')
        reactants = [Chem.MolFromSmiles(x) for x in r_str.split('.')] if r_str else []
        products  = [Chem.MolFromSmiles(x) for x in p_str.split('.')] if p_str else []
        if not reactants or not products or any(m is None for m in reactants + products):
            return None
        rxn = Reactions.ChemicalReaction()
        for m in reactants: rxn.AddReactantTemplate(m)
        for m in products:  rxn.AddProductTemplate(m)
        return rxn
    except Exception:
        return None

def process_file(in_path: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    part = 0
    total_valid = 0
    log.info(f"Processing {in_path} ...")

    # Peek first rows
    try:
        peek = pd.read_csv(in_path, sep='\t', header=0, nrows=3, dtype=str, encoding='latin1', on_bad_lines='skip')
        log.info("Sample:\n" + str(peek))
    except Exception as e:
        log.error(f"Failed to read sample: {e}")
        return

    for chunk in pd.read_csv(in_path, sep='\t', header=0, names=COLS, dtype=str,
                             encoding='latin1', on_bad_lines='skip', chunksize=CHUNKSIZE):
        # Basic cleaning
        chunk = chunk[chunk['ReactionSMILES'].notna()].copy()
        chunk = chunk[chunk['ReactionSMILES'] != 'ReactionSmiles']  # drop header row if present
        # require arrow
        chunk = chunk[chunk['ReactionSMILES'].str.contains('>', regex=False)]
        # drop obvious SMARTS/dummies
        chunk = chunk[~chunk['ReactionSMILES'].str.contains(BAD_TOKENS_RE, na=False)]
        # parse year
        chunk['Year'] = pd.to_numeric(chunk['Year'], errors='coerce')

        if chunk.empty:
            continue

        tqdm.pandas(desc="Parsing reactions")
        rxn_objs = chunk['ReactionSMILES'].progress_apply(safe_rxn_from_string)
        valid = chunk[rxn_objs.notna()].copy()
        if valid.empty:
            continue

        # store canonical text (strings are parquet-friendly; RDKit objects are not)
        def _rxn_smiles(r):
            try: return Reactions.ReactionToSmiles(r)
            except Exception: return None
        def _rxn_smarts(r):
            try: return Reactions.ReactionToSmarts(r)
            except Exception: return None

        valid['ReactionSMILES_canonical'] = rxn_objs[rxn_objs.notna()].map(_rxn_smiles)
        valid['ReactionSMARTS'] = rxn_objs[rxn_objs.notna()].map(_rxn_smarts)
        valid = valid.dropna(subset=['ReactionSMILES_canonical'])

        if valid.empty:
            continue

        out_cols = ['ReactionSMILES','ReactionSMILES_canonical','ReactionSMARTS',
                    'PatentNumber','ParagraphNum','Year','TextMinedYield','CalculatedYield']
        out_path = os.path.join(out_dir, f"{os.path.basename(in_path)}.part{part}.parquet")
        valid[out_cols].to_parquet(out_path, index=False)
        total_valid += len(valid)
        part += 1
        log.info(f"Wrote {len(valid)} rows -> {out_path}")

    log.info(f"Done: {in_path}. Valid reactions: {total_valid}")

def main():
    input_files = [
        "../../../5104873/1976_Sep2016_USPTOgrants_smiles.rsmi",
        "../../../5104873/2001_Sep2016_USPTOapplications_smiles.rsmi",
    ]
    out_dir = "./uspto_parquet"
    for f in input_files:
        if os.path.exists(f):
            process_file(f, out_dir)
        else:
            log.error(f"File not found: {f}")

if __name__ == "__main__":
    main()
