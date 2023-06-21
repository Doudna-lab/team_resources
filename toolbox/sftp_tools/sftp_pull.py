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


def download_directory(sftp, remote_dir, local_dir):
    # Recursively download files from the remote directory to the local directory
    for item in sftp.listdir_attr(remote_dir):
        remote_path = remote_dir + '/' + item.filename
        local_path = local_dir + '/' + item.filename

        if item.st_mode & 0o40000:
            # Directory
            if not os.path.exists(local_path):
                os.makedirs(local_path)
            download_directory(sftp, remote_path, local_path)
        else:
            # File
            sftp.get(remote_path, localpath=local_path, callback=transfer_callback)
