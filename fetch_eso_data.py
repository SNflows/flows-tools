import asyncio
from typing import Collection, Optional
import numpy as np
import pandas as pd
from requests.exceptions import HTTPError
from tendrils import api

eso_sites = [2, 12, 15]
eso_site_names = {2: "VLT-HAWKI", 12: "NTT-SOFI", 15: "NTT-EFOSC"}
usecols = ["fileid_img", "time", "mag", "site", "photfilter"]
tgtcols = ["target_name", "redshift", "ra", "decl"]


def get_targets(sub_targets: Optional[list[int]] = None) -> pd.DataFrame:
    targets = api.get_targets()
    targets = pd.DataFrame(targets)
    targets = targets[targets["target_status"] == "target"]
    if sub_targets is not None:
        targets = targets[targets["targetid"].isin(sub_targets)]
    return targets.set_index("targetid")


def get_target_data(tgtid: int, targets: pd.DataFrame) -> pd.Series:
    return targets.loc[(tgtid, tgtcols)]


def get_sites(fids: Collection[int]) -> pd.DataFrame:
    keys = ["obstime", "site", "photfilter"]
    r = {fid: {k: dat[k] for k in keys} for fid in fids if (dat := api.get_datafile(fid))["site"] in eso_sites}
    df = pd.DataFrame(r).T
    df.index.name = "fileid_img"
    df.rename({"obstime": "time"}, axis=1, inplace=True)
    df["mag"] = np.nan
    return df


async def get_target_lc(tgtid: int, targets: pd.DataFrame) -> pd.DataFrame:
    lc = api.get_lightcurve(tgtid)
    lc["time"] = lc["time"].mjd  # type: ignore
    lc = lc.to_pandas()
    lc["mag"] = lc["mag_sub"].fillna(lc["mag_raw"])
    lc = lc[usecols]
    return lc.set_index("fileid_img")


async def get_missing_eso_lcs(lc: pd.DataFrame, tgtid: int, targets: pd.DataFrame) -> pd.DataFrame:
    # add data with failed/missing photometry
    dfs = api.get_datafiles(tgtid, filt="all")
    missing_fids = set(dfs) - set(lc.index.values)
    lc = pd.concat((lc, get_sites(missing_fids)))
    lc.sort_index(inplace=True)
    return lc[lc["site"].isin(eso_sites)]


async def get_eso_lc(tgtid: int, targets: pd.DataFrame) -> pd.DataFrame:
    try:
        lc = await get_target_lc(tgtid, targets)
    except HTTPError:
        lc = pd.DataFrame(index=[], columns=usecols)
        lc.set_index("fileid_img", inplace=True)
    lc = await get_missing_eso_lcs(lc, tgtid, targets)
    if len(lc) == 0:
        return pd.DataFrame()
    lc[tgtcols] = get_target_data(tgtid, targets)
    return lc


async def main():
    targets = get_targets()
    coroutines = [get_eso_lc(tgtid, targets) for tgtid in targets.index]
    eso_lcs_with_excp = await asyncio.gather(*coroutines, return_exceptions=True)
    eso_lcs = pd.concat([lc for lc in eso_lcs_with_excp if isinstance(lc, pd.DataFrame)])
    eso_lcs["site"] = eso_lcs.site.apply(lambda x: eso_site_names[x])
    eso_lcs.to_json("eso_lcs.json")
    eso_lcs.to_csv("eso_lcs.csv")
    print([lc for lc in eso_lcs_with_excp if not isinstance(lc, pd.DataFrame)])


if __name__ == "__main__":
    asyncio.run(main())
