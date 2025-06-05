import logging
import pandas as pd
from typing import List, Optional, Any
from dataclasses import dataclass, field
import solara

# Set up logging
logger = logging.getLogger(__name__)

@dataclass
class AppState:
    """Reactive state for the application."""
    genes: List[str] = field(default_factory=lambda: ["Azin2", "Tpi1"])
    df: Optional[pd.DataFrame] = field(default_factory=lambda: None)
    adata: Optional[Any] = field(default_factory=lambda: None)
    loading_error: Optional[str] = field(default_factory=lambda: None)
    is_loading: bool = field(default_factory=lambda: False)
    organism: str = field(default_factory=lambda: "Mouse")
    clustering_ids: List[str] = field(default_factory=list)

# Reactive state variables
solara_state = solara.reactive(AppState())
cluster_id = solara.reactive("RNA_snn_res.0.1")  # Standalone reactive for cluster_id
progress = solara.reactive(0.0)  # Standalone reactive for progress

ORGANISMS = ["Mouse", "Human"]

