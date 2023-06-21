import pysftp
import os


def transfer_callback(event, chunk):
    # Callback function to track the transfer progress
    if event == 'put':
        # File upload event
        print(f"Uploaded {chunk} bytes")
    elif event == 'get':
        # File download event
        print(f"Downloaded {chunk} bytes")


def download_directory(sftp_instance, remote_dir, local_dir):
    # Recursively download files from the remote directory to the local directory
    for item in sftp_instance.listdir_attr(remote_dir):
        remote_path = remote_dir + '/' + item.filename
        local_path = local_dir + '/' + item.filename

        if item.st_mode & 0o40000:
            # Directory
            if not os.path.exists(local_path):
                os.makedirs(local_path)
            download_directory(sftp_instance, remote_path, local_path)
        else:
            # File
            sftp_instance.get(remote_path, localpath=local_path, callback=transfer_callback)


# Create an instance of the `pysftp.Connection` class
with pysftp.Connection('sftp.genewiz.com',
                       username='isylvain_berkeley',
                       password='aS06jJvl2oMG0MvKIBtx') as sftp:
    # Change to the desired remote directory
    sftp.chdir('30-857011403/00_fastq/')

    # Start the recursive download of files
    download_directory(sftp, '.', '/groups/doudna/team_resources/azenta_temp/30-857011403')
