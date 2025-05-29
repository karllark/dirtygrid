import os
import glob
import shutil
from astropy.io import fits

if __name__ == '__main__':

    filepath = "/astro/dust_kg3/klaw/cloudy2/nasa_fits"
    outpath = "/astro/dust_kg3/kgordon/DirtyGrid"

    alldirs = glob.glob(f"{filepath}/??") + glob.glob(f"{filepath}/?")
    mdirs = sorted(alldirs)

    for cmdir in mdirs:
        print(cmdir)
        psdirs = glob.glob(f"{cmdir}/*")
        for csdir in psdirs:
            pfiles = glob.glob(f"{csdir}/*.fits")
            for cfile in pfiles:
                newname = f"s{cfile.split("-s")[1]}"
                newname = cfile[cfile.find("-s"):]

                # info from filename
                ggeom = (cfile.split("-")[-1]).split("_")[0]
                if (ggeom == "homo") or (ggeom == "clumpy"):
                    ggeom = (cfile.split("-")[-1]).split("_")[1]

                # info from header
                header = fits.getheader(cfile, ext=1)
                gtypes = [s for s in header["COMMENT"] if "/Models/" in s]
                gtype = ((gtypes[0].split("/"))[-1]).lower()

                ffactors = [s for s in header["COMMENT"] if "density_ratio=" in s]
                ffactor = (ffactors[0].split("="))[-1]
                if "0.01" in ffactors[0]:
                    lgeom = "clumpy"
                    newname = newname.replace(f"-{ggeom}", "")
                else:
                    lgeom = "homo"
                    newname = newname.replace(f"-{lgeom}_{ggeom}", "")

                newname = f"{gtype}_{ggeom}_{lgeom}{newname}"

                newpath = f"{outpath}/{gtype}/{ggeom}/{lgeom}"
                print(f"{newpath}/{newname}", cfile)
                if not os.path.exists(newpath):
                    os.makedirs(newpath)

                shutil.copyfile(cfile, f"{newpath}/{newname}")