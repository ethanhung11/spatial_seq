import pandas as pd
import liana as li

from typing import List, Literal
import requests
import io
import re
from tqdm import tqdm


def get_cellphonedbv5_resource():
    # see https://github.com/saezlab/pypath/issues/269#issuecomment-2538243027
    resource = requests.get(
        "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv"
    ).content
    resource = io.StringIO(resource.decode("utf-8"))
    resource = pd.read_csv(resource, sep=",")
    # keep only PPIs
    resource = resource[resource["is_ppi"]][["interactors"]]
    # replace + with _
    resource["interactors"] = resource["interactors"].apply(
        lambda x: x.replace("+", "_")
    )
    # if interactors contains two '-' replace the first one with '&
    resource["interactors"] = resource["interactors"].apply(
        lambda x: x.replace("-", "&", 1) if x.count("-") == 2 else x
    )
    # split by - and expand
    resource = resource["interactors"].str.split("-", expand=True)
    # replace & with - in the first column
    resource[0] = resource[0].apply(lambda x: x.replace("&", "-"))
    resource.columns = ["ligand", "receptor"]
    return resource


def human2mouse(human, map_df=None):
    if map_df is None:
        map_df = li.rs.get_hcop_orthologs(
            url="https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz",
            columns=["human_symbol", "mouse_symbol"],
            # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
            min_evidence=3,
        )

        # rename the columns to source and target, respectively for the original organism and the target organism
        map_df = map_df.rename(
            columns={"human_symbol": "source", "mouse_symbol": "target"}
        )

    # Translate Resource
    mouse = li.rs.translate_resource(
        human,
        map_df=map_df,
        columns=["ligand", "receptor"],
        replace=True,
        # Here, we will be harsher and only keep mappings that don't map to more than 1 mouse gene
        one_to_many=1,
    )
    return mouse


def liana_mouse_resource(
    resource_name: (
        Literal[
            "baccin2019",
            "cellcall",
            "cellchatdb",
            "cellinker",
            "cellphonedb",
            "cellphondbv5",
            "celltalkdb",
            "connectomedb2020",
            "consensus",
            "embrace",
            "guide2pharma",
            "hpmr",
            "icellnet",
            "italk",
            "kirouac2010",
            "lrdb",
            "mouseconsensus",
            "ramilowski2015",
        ]
        | List[str]
    ),
):
    if type(resource_name) is str:
        resource_name = [resource_name]

    resources = dict()
    map_df = li.rs.get_hcop_orthologs(
        url="https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz",
        columns=["human_symbol", "mouse_symbol"],
        # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
        min_evidence=3,
    )
    map_df = map_df.rename(columns={"human_symbol": "source", "mouse_symbol": "target"})

    for name in (pbar := tqdm(resource_name)):
        pbar.set_description(f"Processing {name}")
        if name == "cellphonedbv5":
            res = get_cellphonedbv5_resource()
        else:
            res = li.rs.select_resource(name)

        isMouse = bool(re.search(r"[a-z]", res.iloc[0, 0]))
        if isMouse is True:
            print(
                f"found lowercase letters in first cell '{res.iloc[0,0]}' of '{name}', not converting!"
            )
        else:
            res = human2mouse(res, map_df)

        resources[name] = res

    if len(resources) > 1:
        dfs = []
        for name, df in resources.items():
            temp = df.drop_duplicates().copy()
            temp["db_sources"] = name
            dfs.append(temp)

        # Concat all dataframes
        merged = pd.concat(dfs, ignore_index=True)
        cols = [c for c in merged.columns if c != "db_sources"]

        # Track duplicate sources
        for _, group in merged.groupby(cols):
            if len(group) > 1:
                sources = ", ".join(sorted(group["db_sources"]))
                merged.loc[group.index, "db_sources"] = sources

        # Remove duplicate rows, keep first
        resources = merged.drop_duplicates(subset=cols, keep="first").reset_index(
            drop=True
        )

    else:
        resources = resources[name]

    return resources
