#!/usr/bin/env python
"""
Example of how to make a custom parser.
"""

import vcf_parser_dev as vp


class TestParse(vp.Parser):
    pass


vp.main()

