# phylocorrelate

## Dependencies

- `BLAST+`
- `R` version 3.6.x or later

## R package dependencies

### CRAN

`shiny`, `shinyWidgets`, `DT`, `shinycssloaders`, `igraph`, `plotly`, `visNetwork`, `shinythemes`, `ggplot2`, `reshape2`, `magrittr`, `promises`, `future`, `readr`, `gmailr`, `ps`

### BioConductor (version 10 or later)

`topGO`, `S4Vectors`

## Gettings emails from the community message function

Right now the apps send emails automatically whenever someone submits a message using the community feature. To do this, it uses local GMail keys. For obvious reasons, these are not included in the GitHub repo.

## Configuration

Inside the `app` folder is a configuration file, `PhyloCorrConfig.txt`. Editing this can change some of the behaviours of the app. Remember to change the file paths to match your system.

- `URL = "https://phylocorrelate.uwaterloo.ca/app"`: server app URL
- `IndexURL = "https://phylocorrelate.uwaterloo.ca"`: server index URL
- `UseBlastp = FALSE`: `TRUE`/`FALSE`; enable/disable BLAST functionality
- `BlastpDB = "/home/$USER/phylocorrelate/blastdb/AllProteomes.faa"`: Location of BLAST database
- `UseGmailBlastp = FALSE`: `TRUE`/`FALSE`; enable/disable the server sending email for BLAST jobs via GMail
- `ServerEmail = "phylocorrelate@gmail.com"`: email to be used by the server to send emails
- `AdminEmail = "phylocorrelate@gmail.com"`: admin email; the server will send emails to this address in case of errors during BLAST analyses
- `GmailCredentials = "/home/$USER/phylocorrelate/app/blastp/credentials.json"`: GMail credentials; look up the [`gmailr`](https://github.com/r-lib/gmailr) package for instructions on how to set this up
- `GmailCache = "/home/$USER/phylocorrelate/app/blastp/.secret"`: GMail cache; look up the [`gmailr`](https://github.com/r-lib/gmailr) R package for instructions on how to set this up
- `BlastpThreads = 2`: number of threads to use when running BLAST

## Generating the necessary data files and BLAST database

**Warning: not recommended unless you have a supercomputer or are willing to wait weeks for the analyses to complete.**

Instructions on how to download and analyse the PFAM/TIGRFAM/KO tables go here.

## Running the app locally

```r
shiny::runApp("app")
```

## Running the app with Shiny Server open source edition

Note: these instructions were written with Ubuntu 16.04.6 LTS in mind. They also assume you've created the log folders.

### No SSL

Edit your `/etc/shiny-server/shiny-server.conf` (replace `$USER` with your username):

```
run_as shiny;

sanitize_errors false;
preserve_logs true;

server {
  listen 80;
  location / {
      run_as $USER;
      directory_index on;
      app_idle_timeout 0;
      site_dir /home/$USER/phylocorrelate/index/;
      log_dir /home/$USER/phylocorrelate-logs/index/;
    }
  location /app {
      run_as $USER;
      app_idle_timeout 0;
      app_dir /home/$USER/phylocorrelate/app/;
      log_dir /home/$USER/phylocorrelate-logs/app/;
    }
}
```

If you ever make changes to the app, restart the server:

```
sudo systemctl restart shiny-server
```

### With SSL

If you've installed Shiny Server already, stop it:

```
sudo systemctl stop shiny-server
```

Install `nginx` and stop it:

```
sudo apt-get install nginx
sudo systemctl stop nginx
```

Now edit their configs:

```
# /etc/nginx/sites-available/default
server {
  listen 80;
  return 301 https://$host$request_uri;
}
server {
  listen 443 ssl;
  server_name phylocorrelate.uwaterloo.ca;
  ssl_certificate /etc/letsencrypt/live/phylocorrelate.uwaterloo.ca/fullchain.pem;
  ssl_certificate_key /etc/letsencrypt/live/phylocorrelate.uwaterloo.ca/privkey.pem;
  location / {
    proxy_pass http://127.0.0.1:3838;
  }
}
```

Replace `$USER` with the server username:

```
# /etc/shiny-server/shiny-server.conf
run_as shiny;
sanitize_errors false;
preserve_logs true;
server {
  listen 3838 127.0.0.1;
  location / {
    run_as $USER;
    directory_index on;
    app_idle_timeout 0;
    site_dir /home/$USER/phylocorrelate/index/;
    log_dir /home/$USER/phylocorrelate-logs/index/;
  }
  location /app-A {
    run_as $USER;
    app_idle_timeout 0;
    app_dir /home/$USER/phylocorrelate/app/;
    log_dir /home/$USER/phylocorrelate-logs/app/;
  }
}
```

Install `certbot` (based on these [instructions](https://geekflare.com/setup-nginx-with-lets-encrypt-cert/)):

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:certbot/certbot
sudo apt-get update
sudo apt-get install python-certbot-nginx
```

Run `certbot`:

```
sudo certbot --nginx
```

It will ask for a domain name. For example, `diffbase.uwaterloo.ca`. Once it's done it will tell you when the certificate will expire. When the time comes, renew it with:

```
sudo certbot renew
```

The `certbot` program will try and automatically modify your `nginx` config file. Make sure it did it correctly.

Start `nginx` and `shiny-server`:

```
sudo systemctl start nginx
sudo systemctl start shiny-server
```
