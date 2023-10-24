## Introduction
- This tutorial page covers the necessary steps to conduct data intake from Azenta and its integration with BaseSpace 
- The placeholder paths and credentials used in this tutorial are:
  - `<Azenta SFTP login>`: `isylvain_berkeley@sftp.genewiz.com`
  - `<remote sftp path>`: `30-857011403/00_fastq/`
  - `<local path to store download>`: `/groups/doudna/team_resources/azenta_temp/30-857011403`
  - `<BaseSpace project label>`: `IGI_Germline_WGS_Assay`
  - `<BaseSpace project id>`: `392083693`

## BaseSpace Setup
### Local machine setup
`bs auth --force`

### Copy URL to the browser while logged on BaseSpace 'IGI Clinical Lab' account
### Double-check if all the projects under the account show up 
`bs list projects`

## SFTP Download
### SFTP Login
`sftp isylvain_berkeley@sftp.genewiz.com`

## Check the local working directory (where the data will be downloaded)
`lpwd`

### The outcome should be:
### Local working directory: /groups/clinical/projects/Azenta_Data/30-857011403
### If not: Setup the local directory for download
`lcd /groups/clinical/projects/Azenta_Data/30-857011403`    

### Go to the directory where the content is located
`cd 30-857011403/00_fastq/`

### Proceed with the download
`mget *`

### Some downloads can be too big and should be downloaded on the background to avoid network interruptions
### The script /groups/doudna/team_resources/toolbox/sftp_tools/sftp_pull.py can work around that
### Further instructions will be available soon


## Adjust FASTQ name convention 
### BaseSpace requires FASTQ filenames to adhere to a specific naming convention
### The script fastq_refactor.py gathers all filenames within a directory and apply
### the BaseSpace designated rules
`cd /groups/doudna/team_resources/toolbox/baseSpace_tools`
`python3 fastq_refactor.py /groups/clinical/projects/Azenta_Data/30-857011403`

## Upload to BaseSpace
### Create project in the BaseSpace account
`bs create project -n IGI_Germline_WGS_Assay`

### Go to the local directory where the content is located
`cd /groups/clinical/projects/Azenta_Data/30-857011403`

### Upload the dataset
`nohup bs upload dataset --project=392083693 /groups/clinical/projects/Azenta_Data/30-857011403/* &`
