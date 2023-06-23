
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
### Local working directory: /groups/doudna/team_resources/azenta_temp/30-857011403
### If not: Setup the local directory for download
`lcd /groups/doudna/team_resources/azenta_temp/30-857011403`    

### Go to the directory where the content is located
`cd 30-857011403/00_fastq/`

### Proceed with the download
`mget *`

### Some downloads can be too big and should be downloaded on the background to avoid network interruptions
### The script /groups/doudna/team_resources/toolbox/sftp_tools/sftp_pull.py can work around that
### Further instructions will be available soon


## Upload to BaseSpace
### Create project in the BaseSpace account
`bs create project -n 30-857011403 -d "project description"`

### Go to the directory where the content is located
`cd /groups/doudna/team_resources/azenta_temp/30-857011403`

### Upload the dataset
`bs upload dataset --project=30-857011403 *`

