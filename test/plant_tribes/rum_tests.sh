#!/bin/sh

# Usage: test.sh [--all|--assembly_post_processor|--gene_family_classifier|
#                 --gene_family_integrator|--gene_family_aligner|
#                 --gene_family_phylogeny_builder|--kaks_analysis|
#                 --ks_distribution]
# NOTE: The kaks_analysis tool test hangs on debruijn due to a dependency
# permissions issue when the tool is executed within this test framework.
# The tool does not hang when run from within Galaxy.  This issue is why
# the --al option below does not include the --kaks_analysis flag.

for arg in "$@"; do
    case "$arg" in
        --all)
            /bin/bash ./test.sh --assembly_post_processor --gene_family_classifier --gene_family_integrator --gene_family_aligner --gene_family_phylogeny_builder --ks_distribution
            ;;

        --assembly_post_processor)
            /bin/bash ./test.sh --assembly_post_processor
            ;;

        --gene_family_classifier)
            /bin/bash ./test.sh --gene_family_classifier
            ;;

        --gene_family_integrator)
            /bin/bash ./test.sh --gene_family_integrator
            ;;

        --gene_family_aligner)
            /bin/bash ./test.sh --gene_family_aligner
            ;;

        --gene_family_phylogeny_builder)
            /bin/bash ./test.sh --gene_family_phylogeny_builder
            ;;

        --kaks_analysis)
            /bin/bash ./test.sh --kaks_analysis
            ;;

        --ks_distribution)
            /bin/bash ./test.sh --ks_distribution
            ;;
    esac
done

