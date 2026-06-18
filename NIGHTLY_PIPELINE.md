# Nightly Pipeline Setup

## Overview
Automated nightly data downloads and processing pipeline with email alerting.

## Scripts Location
All scripts located in: `/home/cfranken/scripts/`

## Components

### 1. Data Synchronization
**Script:** `sat_sync.py`
- Downloads IDP, Met, MODIS data from DAAC
- Configuration: `sat_sync.toml`
- Cron: 2:30 AM daily (`cron_sat_sync_and_merge.sh`)

### 2. Data Merging
**Script:** `merge_imap_aband_netcdf.py`
- Merges IDP + Met into single NetCDF per orbit
- Output: `/kiwi-data/Data/satellite/OCO2/IMAP-SIF/merged/`
- Runs after sat_sync completes

### 3. L1b SNR Subsetting
**Script:** `download_l1b_snr_subset.py`
- Downloads SNR coefficients + radiance subset from L1b files
- OPeNDAP subsetting: `snr_coef`, `MaxMS`, `radiance_o2/strong/weak[5:1016:100]`
- Output: `/kiwi-data/Data/satellite/OCO2/L1b_SNR/<year>/<doy>/`
- Cron: 3:00 AM daily (`cron_download_l1b_snr.sh`)

**Extraction Tools:**
- `extract_snr_radiance_subset.py` — Extract SNR from local or OPeNDAP L1b
- `extract_snr_every10th.py` — Historical timeseries (earthaccess method)
- `extract_snr_direct_opendap.py` — Direct OPeNDAP access (requires auth)

### 4. Health Monitoring
**Script:** `check_nightly_logs.py`
- Monitors all nightly run logs
- Detects errors, warnings, failed downloads
- Sends email alert if issues found
- Runs automatically after main pipelines (~3:30 AM)

## Cron Schedule
```
30 2 * * * /home/cfranken/scripts/cron_sat_sync_and_merge.sh    # Sync + merge
0  3 * * * /home/cfranken/scripts/cron_download_l1b_snr.sh      # L1b SNR download
→ ~3:30 AM: Both scripts + log checks complete, email sent if needed
```

## Log Files
Located in: `/home/cfranken/scripts/logs/`

- `sat_sync_YYYYMMDD.log` — IDP/Met/MODIS downloads
- `merge_YYYYMMDD.log` — NetCDF merging
- `l1b_snr_download_YYYYMMDD.log` — L1b SNR downloads
- `cron_YYYYMMDD.log` — Overall orchestration

## Output Data

### Merged IDP + Met
```
/kiwi-data/Data/satellite/OCO2/IMAP-SIF/merged/
oco2_L2MERGED{mode}_{orbit}_{date}_{time}.nc4
```
Variables: 47 (37 IDP + 10 Met)

### L1b SNR Subsets
```
/kiwi-data/Data/satellite/OCO2/L1b_SNR/<year>/<doy>/
oco2_L1bSc{ND|GL|XS}_*_snr.nc
```
Variables: snr_coef, MaxMS, radiance_o2, radiance_strong_co2, radiance_weak_co2

### Extracted Timeseries (from `extract_*` scripts)
```
/tmp/snr_timeseries.csv
/tmp/radiance_timeseries.csv
```

## Testing

### Test log checker
```bash
python3 /home/cfranken/scripts/check_nightly_logs.py \
  --email cfranken@caltech.edu \
  --date YYYYMMDD \
  --send-always
```

### Test L1b SNR extraction (local files)
```bash
python3 /home/cfranken/scripts/extract_snr_radiance_subset.py --local --limit 1
```

### Test sat_sync manually
```bash
python3 /home/cfranken/scripts/sat_sync.py config \
  /home/cfranken/scripts/sat_sync.toml
```

## Troubleshooting

### Email not sending
- Check SMTP server: `echo | telnet localhost 25`
- Verify email address in scripts (default: cfranken@caltech.edu)
- Check `/home/cfranken/scripts/logs/` for error messages

### L1b SNR downloads failing
- OPeNDAP subsetting issues (netCDF4 constraints)
- Fallback: Use local extraction on downloaded L1b files
- Check Earthdata credentials in `~/.netrc`

### No log files found
- Verify cron job is scheduled: `crontab -l`
- Check log directory exists: `/home/cfranken/scripts/logs/`
- Run scripts manually to generate logs

## Files Modified
- `/home/cfranken/scripts/cron_sat_sync_and_merge.sh` — Added log checker
- `/home/cfranken/scripts/cron_download_l1b_snr.sh` — Added log checker

## Files Created
- `check_nightly_logs.py` — Log monitoring + email alerts
- `download_l1b_snr_subset.py` — L1b SNR subset downloader
- `extract_snr_radiance_subset.py` — SNR extraction from L1b (local/OPeNDAP)
- `extract_snr_every10th.py` — Historical SNR timeseries (earthaccess)
- `extract_snr_direct_opendap.py` — Direct OPeNDAP SNR access

All scripts are production-ready and have been tested.
