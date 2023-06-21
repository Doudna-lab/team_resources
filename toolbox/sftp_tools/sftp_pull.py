import pysftp


def transfer_callback(event, chunk):
    # Callback function to track the transfer progress
    if event == 'put':
        # File upload event
        print(f"Uploaded {chunk} bytes")
    elif event == 'get':
        # File download event
        print(f"Downloaded {chunk} bytes")


# Create an instance of the `pysftp.Connection` class
with pysftp.Connection('sftp.genewiz.com',
                       username='isylvain_berkeley',
                       password='aS06jJvl2oMG0MvKIBtx') as sftp:

    # Change to the desired remote directory
    sftp.chdir('30-857011403/00_fastq/')

    # Start the background file transfer
    sftp.get('30-857011403/00_fastq/', localpath='/groups/doudna/team_resources/azenta_temp/30-857011403', callback=transfer_callback)

# The file transfer is running in the background after the `with` block is exited
