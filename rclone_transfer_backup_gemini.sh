#!/bin/bash

rclone copy --progress --transfers 12 --checkers 12 --copy-links ../landscape_of_RNA_isoform_expression/ gemini3:/mnt/gemini3-4/PROJECTS/landscape_of_RNA_isoform_expression
