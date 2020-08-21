import sys
import requests
import time
import threading
import _thread as thread

def get_csrf(REQ_URL):
    """Get a csrf token from the request url to authenticate further requests
    Parameters
    ----------
    REQ_URL string
        The URL that you want to make a request against after getting the token

    Returns
    -------
    csrftoken
        csrf token to use for further requests against the same URL
    """
    client = requests.session()
    # Retrieve the CSRF token first
    client.get(REQ_URL)  # sets cookie
    if 'csrftoken' in client.cookies:
        # Django 1.6 and up
        csrftoken = client.cookies['csrftoken']
    else:
        # older versions
        csrftoken = client.cookies['csrf']
    return csrftoken


def update_cset(REQ_URL, target_name, submit_choice, upload_key, update_set, sdf_path, pdb_zip_path=None):
    """Send data to <root_url>/viewer/upload_cset/ to update an existing computed set. The existing
    data will be removed, and the new data added. Any molecules not re-uploaded will disappear.

    Parameters
    ----------
    REQ_URL str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/)
    target_name str
        the name of the target in Fragalysis that the computed set is for
    submit_choice int
        0 for validate, 1 for upload
    upload_key str
        upload key, not currently turned on, so can be any value, but not blank or null
    update_set str
        the name of the computed set you want to update
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file)
    sdf_path str
        path to the sdf file to upload
    pdb_zip_path str
        path to the zip file of pdb's to upload (optional)

    Returns
    -------
    taskurl str
        the URL to check for the status of the upload
    """
    print(f'Submitting files to update {update_set}...')

    csrftoken = get_csrf(REQ_URL)

    payload = {'target_name': target_name,
               'submit_choice': submit_choice,
               'upload_key': upload_key,
               'update_set': update_set}

    files = [
        ('sdf_file', open(sdf_path, 'rb')),

    ]

    if pdb_zip_path:
        files.append(('pdb_zip', open(pdb_zip_path, 'rb')))

    headers = {'X-CSRFToken': csrftoken,
               'Cookie': f'csrftoken={csrftoken}'}

    response = requests.request("POST", REQ_URL, headers=headers, data=payload, files=files)

    lines = response.text.split('\n')
    taskurl = None
    for l in lines:
        if 'taskUrl = "/viewer/upload_task/' in l:
            taskid = l.split('/')[-2]
            print(f'upload task id: {taskid}')
            taskurl = f'{REQ_URL.replace("/viewer/upload_cset/", "/viewer/upload_task/")}{taskid}'

            if taskurl:
                break
    if not taskurl:
        raise Exception('Something went wrong with the upload request! \
                        Please try again or email rachael.skyner@diamond.ac.uk for help.')

    return taskurl


def quit_function(fn_name):
    """Quit a function and return an error
    Parameters
    ----------
    fn_name
        name of function to apply to
    """
    # print to stderr, unbuffered in Python 2.
    print('{0} took too long. The task has probably not worked, but is left in a PENDING state. \
    Please try again or email rachael.skyner@diamond.ac.uk for help.'.format(fn_name), file=sys.stderr)
    sys.stderr.flush()  # Python 3 stderr is likely buffered.
    thread.interrupt_main()  # raises KeyboardInterrupt


def exit_after(s):
    '''
    use as decorator to exit process if function takes longer than s seconds

    Parameters
    ----------
    s int
        number of seconds to exit after
    '''

    def outer(fn):
        def inner(*args, **kwargs):
            timer = threading.Timer(s, quit_function, args=[fn.__name__])
            timer.start()
            try:
                result = fn(*args, **kwargs)
            finally:
                timer.cancel()
            return result

        return inner

    return outer


@exit_after(600)
def get_task_response(taskurl):
    """Check a task url to get it's status. Will return SUCCESS or FAILED, or timeout after 10 minutes (600s)
       if the task is still pending

    Parameters
    ----------
    taskurl str
        URL to ping

    Returns
    -------
    status str
        SUCCESS or FAILED
    """
    print('pinging task to check status...')
    requests.request("GET", taskurl)
    complete = False
    while not complete:
        task_response = requests.request("GET", taskurl)
        status = task_response.json()['upload_task_status']
        if status == "SUCCESS":
            complete = True
        if status == "FAILED":
            complete = True
        time.sleep(5)
    return status

# EXAMPLE:
# REQ_URL = 'https://fagalysis.diamond.ac.uk/viewer/upload_cset/'
# taskurl = update_cset(REQ_URL=REQ_URL,
#                       target_name='Mpro',
#                       submit_choice='1',
#                       upload_key='1',
#                       update_set='WT-xCOS3-ThreeHop',
#                       sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28.sdf',
#                       pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs.zip')
# task_response = get_task_response(taskurl)
# print(task_response)