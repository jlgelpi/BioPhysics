#! /usr/bin/python3

"""Example script for building custom command-lines use argparse"""

import argparse
import sys

parser = argparse.ArgumentParser(
                                 prog='ProgName', 
                                 description='Description of the program'
                                 )

# Use add_argument for all needed arguments
# https://docs.python.org/3/library/argparse.html

parser.add_argument(
                    '--option1', 
                    action='store_true', 
                    dest='var1',
                    help='ON/OFF option1'
                    )

parser.add_argument(
                    '--option2', 
                    dest='var2',
                    help='Introduce option2'
                    )

parser.add_argument('required_file',
                    help='Required file for the program')

# Read command line into args

args = parser.parse_args()
       
    
print ("\nSettings\n--------")

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print ("\nSettings, again\n---------------")    

Variable1 = args.var1
Variable2 = args.var2
Variable_Required = args.required_file

print(Variable1, Variable2, Variable_Required)    

