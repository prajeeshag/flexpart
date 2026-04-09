import xarray as xr


def convert_gfs_grib_to_nc(grib_file, nc_file):
    print(f"Opening {grib_file}...")
    ds = xr.open_dataset(
        grib_file, engine="cfgrib", filter_by_keys={"typeOfLevel": "isobaricInPa"}
    )

    print(ds)


if __name__ == "__main__":
    print("Converting GFS GRIB to NC...")
    grib_file = "inputs/gfs.t00z.pgrb2.0p25.f000"
    nc_file = "path/to/gfs.nc"
    convert_gfs_grib_to_nc(grib_file, nc_file)
