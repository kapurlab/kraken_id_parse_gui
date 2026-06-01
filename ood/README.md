# OOD App Deployment

Source of truth for the OOD apps lives here. Changes must be manually copied to the
deployed location (requires sudo). OOD reloads templates on next session launch.

## App structure

```
ood/apps/
  kraken_id_parse_gui/        ← production app
  kraken_id_parse_gui_dev/    ← dev branch-picker app
```

## Deploying (run from repo root)

```bash
# Production app
sudo rsync -av ood/apps/kraken_id_parse_gui/ \
    /var/www/ood/apps/sys/kraken_id_parse_gui/

# Dev app
sudo rsync -av ood/apps/kraken_id_parse_gui_dev/ \
    /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
```

Or copy file by file if rsync isn't available:

```bash
sudo mkdir -p /var/www/ood/apps/sys/kraken_id_parse_gui/template
sudo cp ood/apps/kraken_id_parse_gui/manifest.yml     /var/www/ood/apps/sys/kraken_id_parse_gui/
sudo cp ood/apps/kraken_id_parse_gui/form.yml         /var/www/ood/apps/sys/kraken_id_parse_gui/
sudo cp ood/apps/kraken_id_parse_gui/submit.yml.erb   /var/www/ood/apps/sys/kraken_id_parse_gui/
sudo cp ood/apps/kraken_id_parse_gui/view.html.erb    /var/www/ood/apps/sys/kraken_id_parse_gui/
sudo cp ood/apps/kraken_id_parse_gui/template/before.sh     /var/www/ood/apps/sys/kraken_id_parse_gui/template/
sudo cp ood/apps/kraken_id_parse_gui/template/script.sh.erb /var/www/ood/apps/sys/kraken_id_parse_gui/template/

sudo mkdir -p /var/www/ood/apps/sys/kraken_id_parse_gui_dev/template
sudo cp ood/apps/kraken_id_parse_gui_dev/manifest.yml      /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
sudo cp ood/apps/kraken_id_parse_gui_dev/form.yml          /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
sudo cp ood/apps/kraken_id_parse_gui_dev/submit.yml.erb    /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
sudo cp ood/apps/kraken_id_parse_gui_dev/view.html.erb     /var/www/ood/apps/sys/kraken_id_parse_gui_dev/
sudo cp ood/apps/kraken_id_parse_gui_dev/template/before.sh.erb  /var/www/ood/apps/sys/kraken_id_parse_gui_dev/template/
sudo cp ood/apps/kraken_id_parse_gui_dev/template/script.sh.erb  /var/www/ood/apps/sys/kraken_id_parse_gui_dev/template/
```

## Conda environment setup

Install the shared conda env (do once, requires conda/mamba):

```bash
cd /srv/kapurlab/tools/kraken_id_parse_gui
conda env create -p ./env -f conda_setup/environment.yml

# Install additional bioinformatics tools
conda install -p ./env -c conda-forge -c bioconda \
    blast bwa spades seqkit bracken pigz parallel
```

The OOD script.sh.erb looks for the env at `/srv/kapurlab/tools/kraken_id_parse_gui/env`.

## Frontend build

After any frontend edit:

```bash
cd /srv/kapurlab/tools/kraken_id_parse_gui/frontend
npm install   # first time only
npm run build
```

The compiled output in `frontend/dist/` is served by FastAPI.

## Dashboard — activate the Kraken card

Update `/etc/ood/config/wgs_pipelines.yml` to set the Kraken entry to `status: available`
and add the `launch_url`:

```yaml
  - id: kraken
    icon: k
    title: Kraken
    description: >-
      Read-level taxonomic classification against custom databases. For
      contamination screening and metagenomic profiling.
    status: available
    runtime: "GUI-based interactive fully contained sessions"
    launch_url: /pun/sys/dashboard/batch_connect/sys/kraken_id_parse_gui/session_contexts/new
```
