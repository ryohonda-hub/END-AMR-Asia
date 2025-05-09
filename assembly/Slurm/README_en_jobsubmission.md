# How to Prepare and Run a Job Submission Script  
for Slurm on NIG Supercomputer, Last updated: 2025-04-19

## Preparing the Log Directory  
Create a directory `~/log/` to save job output logs in your home directory on the supercomputer. This only needs to be done once.
```bash
$ mkdir ~/log/
```
---

## Preparing a Job Submission Script  
You can create and edit a script using a source-code editor on your local PC. (You can also use `emacs` on the supercomputer.)
### Example Script 
An example script to submit the job on the execution script`test.sh`:
```bash
#!/bin/bash
#SBATCH --output=/home/your_username/log/%x_%j.out.log
#SBATCH --error=/home/your_username/log/%x_%j.err.log
#SBATCH -p epyc
#SBATCH -t 3-23:59:59
#SBATCH -N 1-1 
#SBATCH -n 24
#SBATCH --mem=48G

./test.sh
```

### Editing Header Lines (lines starting with `#`)  
#### 1. Change the log output directory  
Update lines 2–3 with your own log directory path under your home directory (e.g., `/home/ryohonda/log/%x_%j.out.log`).

#### 2. Specify job execution options (resource requirements)  
From line 4 onward, request the necessary supercomputer resources using `#SBATCH` options.  
You can omit options you don't need. Request appropriate resources depending on the job you will execute.

##### Common Options  
- `-p` *compute node* (i.e,, `epyc`, `short`, `rome`, or `medium`)
- `-t` *max runtime* (`d-hh:mm:ss`). — If omitted, the default `TIMELIMIT` for the node applies (e.g., 1 hour for the *short* node).  
- `-N` *min-max number of nodes* — Specify only for a parallel job (when you request multiple CPU threads). Typically `-N 1-1`.  
- `-n` *number of CPU threads*  — Must specify `-N` as well when you request multiple threads.  
- `--mem=` *total memory for the job* (e.g., `--mem=32G`)
- `-J` *job name (optional)*

### Specify the Job Execution Script(s) to Be Run  
After the header lines, add the commands or job execution scripts you want to run.  
If the execution script is in the same directory as the submission script, use `./` (meaning "current directory").

You can execute multiple scripts sequentially as follows:
```bash
#!/bin/bash
#SBATCH --output=/home/your_username/log/%x_%j.out.log
#SBATCH --error=/home/your_username/log/%x_%j.err.log

./test1.sh
./test2.sh
./test3.sh
```

### Saving and Transferring the Submission Script  
Save the script with `.sh` as the file extension.
Make sure to save the script with **LF (Line Feed)** line endings.  
(Windows uses "CRLF", which can cause errors—**be sure to convert it**!)  

After saving it on your local PC, upload it to the supercomputer via FTP. 

From then on, all work will be done on the supercomputer.

---

## Permission Settings  
### Check File Permissions  
Navigate to the directory with your script and run:
```bash
$ ls -l
```
If the fourth character is `x` (e.g., `-rwxr--r--`), the file is executable.
By default, regular files have the permission `-rw-r--r--`, which does not include execute permission (`x`). Follow the steps below to change the permission of the script you created.

### Add Execute Permission  
Before running, add execute permission for the script owner (`u`):
You need to do this for **both the job submission script and the job execution script**.
```bash
$ chmod u+x script_name.sh
```
To change permissions for all `.sh` files at once:
```bash
$ chmod u+x *.sh
```
Afterward, confirm that both the submission and execution scripts have execute permission (`x`) using `ls -l`.

## Submitting the Job  
*\*Submit jobs from the interactive node of the supercomputer* (login by `$ ssh a001` from the gateway node).

Submit your job submission script using the `sbatch` command.
Example - Submitting the job submission script `js_test.sh`:
```bash
$ sbatch js_test.sh
```

### Checking Job Status  
Check the status of your job with the following command.  
If the status is `PD` or `R`, your job was submitted successfully:
```bash
$ squeue -u your_username
```
When your job completes (status `R`), it disappears from the list.  
You can check output in the log files afterward.

If the job does **not appear in the list** right after submission, it may have failed.  
Check the log files for errors.

#### Job Status Codes  
- **PD** — *PENDING* (waiting for resources)  
- **R** — *RUNNING*  
- **CD** — *COMPLETED* (finished successfully)  
- **CA** — *CANCELLED* (cancelled by user or admin)  
- **F** — *FAILED* (error occurred)  
- **TO** — *TIMEOUT* (ran over the time limit)

### Canceling a Job  
You can cancel a running or pending job with this command (check job ID via `squeue` as above):
```bash
$ scancel job_id
```

## Checking Output and Errors  
Standard output and errors are saved in the `~/log/` directory.  
After job completion, always check for errors.

You can view logs even while the job is running to monitor progress or spot errors.

Use `ls -lt` in the `log` directory to list log files by most recent:
```bash
$ cd ~/log/
$ ls -lt
```

### Types of Log Files  
Two types of log files are generated (in plain text):
- `.out.log` — standard output from scripts and programs (e.g., job progress)  
- `.err.log` — errors or warnings from the system or programs

In the `ls -lt` output, the number before the date shows the file size.  
If the `.err.log` log is 0 bytes, no errors occurred (you can tell without opening it).

### Viewing Log Files with `more`  
Use the `more` command to view log content without downloading it:
```bash
$ more filename
```
If the content is long, you'll see a reversed text line like `--More--(xx%)`.  
Press the spacebar to scroll, or type `q` to quit.

`more` is useful not only for log files, but for any text file (e.g., scripts or BLAST output).  
For more advanced viewing, you can also use `less`.
