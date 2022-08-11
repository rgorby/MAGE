#!/usr/bin/env python

import os
import sys
import subprocess
import time
import argparse
from argparse import RawTextHelpFormatter

from os.path import expanduser
sys.path.insert(1, "./python-slackclient")
from slack import WebClient
from slack.errors import SlackApiError
import logging
logging.basicConfig(level=logging.DEBUG,filename='pyunitReport.log')

def postWithSlackBot(client,message,logname='',channel='#kaijudev'):
    try:
        if (logname):
            response = client.files_upload(
                file=logname,
                initial_comment=message,
                channels=channel,
                )
        else:
           response = client.chat_postMessage(
                text=message,
                channel=channel
           ) 

        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'


if __name__ == '__main__':
    MainS = """Reports on the results of the python unit test."""

    parser = argparse.ArgumentParser(description=MainS,
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t',action='store_true',default=False,
        help='Enables testing mode')
    parser.add_argument('-l',action='store_true',default=False,
        help='Enables loud mode')
    parser.add_argument('-a','--account',type=str,default='UJHB0015',
        help='qsub account number')

    args = parser.parse_args()
    isTest = args.t 
    beLoud = args.l 
    account = args.account 

    # Get Slack API token
    slack_token = os.environ["SLACK_BOT_TOKEN"]
    logging.debug('%s',slack_token)
    client = WebClient(token=slack_token)

    # Get CWD and set main kaiju folder to "home"
    calledFrom = os.path.dirname(os.path.abspath(__file__))
    os.chdir(calledFrom)
    orig = os.getcwd()
    os.chdir('..')
    home = os.getcwd()

    # Go pytest folder
    work = os.path.join(home,'pytests')
    os.chdir(work)

    #remove old results and run pytest
    logname = 'kaiju-pyunit.txt'
    if os.path.exists(logname):
        with open(logname) as f:
            for line in f:
                pass
            last_line = line
    else:
        if(not isTest and (beLoud)):
            postWithSlackBot(client,
                'Python Unit Tests job did not complete\n\n')
            sys.exit(0)                    
        else:
            print("Python Unit Tests job did not complete")
            sys.exit(0)

    hasFail = False
    hasError = False
    hasPass = False 
    if 'fail' in last_line:
        hasFail = True
    if 'error' in last_line:
        hasError = True
    if 'passed' in last_line:
        hasPass = True

    if (isTest): # if isTest do not use slack
        if hasError:
            print("Python Unit Tests did not run sucessfully")
        if hasFail:
            print("Python Unit Tests Failed")
        if (hasPass and not hasError and not hasFail):
            print("Python Unit Tests Passed")
        if (not hasPass and not hasError and not hasFail):
            print('Unexpected error occured during python unit tests\n\n')
    elif (beLoud): # if beLoud send all messages to Slack    
        if hasError:
            postWithSlackBot(client,'Python Unit Tests did not run succesfully\n\n')            
        if hasFail:
            postWithSlackBot(client,'Python Unit Tests Failed\n\n')            
        if (hasPass and not hasError and not hasFail):
            postWithSlackBot(client,'Python Unit Tests Pass')            
        if (not hasPass and not hasError and not hasFail):
            postWithSlackBot(client,'Unexpected error occured during python unit tests\n\n')            
    else: #only report errors to Slack
        if hasError:
            postWithSlackBot(client,'Python Unit Tests did not run succesfully\n\n')            
        if hasFail:
            postWithSlackBot(client,'Python Unit Tests Failed\n\n')           
        if (not hasPass and not hasError and not hasFail):
            postWithSlackBot(client,'Unexpected error occured during python unit tests\n\n')

