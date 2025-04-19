#!/bin/env python

from ete3 import NCBITaxa

def main():
	ncbi = NCBITaxa()
	ncbi.update_taxonomy_database()

if __name__ == '__main__':
    main()