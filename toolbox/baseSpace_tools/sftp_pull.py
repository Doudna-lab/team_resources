import pysftp
import os
import yaml

# Load sftp_config file
with open('config/basespace_sftp.yaml', "r") as f:
    sftp_config = yaml.safe_load(f)


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
        print(f"Processing item {item}")

        if item.st_mode & 0o40000:
            # Directory
            if not os.path.exists(local_path):
                os.makedirs(local_path)
            download_directory(sftp_instance, remote_path, local_path)
        else:
            # File
            sftp_instance.get(remote_path, localpath=local_path, callback=transfer_callback)


try:
    # Create an instance of the `pysftp.Connection` class
    with pysftp.Connection(sftp_config["sftp_address"], sftp_config["username"], sftp_config["password"]) as sftp:
        # Change to the desired remote directory
        sftp.chdir(sftp_config["sftp_target_directory"])

        # Start the recursive download of files
        download_directory(sftp, '.', sftp_config["local_download_directory"])

except Exception as e:
    print(f"An error occurred: {e}")
