#!/bin/sh

# Usage: test.sh [--assembly_post_processor|--gene_family_classifier|
#                 --gene_family_integrator|--gene_family_aligner|
#                 --gene_family_phylogeny_builder|--kaks_analysis|
#                 --ks_distribution]

# Make sure we're running within the appropriate virtual environment.
. /scratch/biotools/galaxy/test/.venv/bin/activate

# Set up the environment
CONDA_CHANNELS='iuc,bioconda,conda-forge,defaults'
CONDA_EXEC='/scratch/biotools/galaxy/_conda/bin/conda'
CONDA_PREFIX='/scratch/biotools/galaxy/_conda'
GALAXY_ROOT_PRISTINE='/scratch/biotools/galaxy/test/galaxy_pristine'
GALAXY_ROOT_TEST='/scratch/biotools/galaxy/galaxy'
LOG_DIR='/scratch/biotools/galaxy/test/log'
LOG_FILE=$LOG_DIR/test.log
PLANTTRIBES_TEST_ROOT='/scratch/biotools/galaxy/test/galaxy_tools/tools/phylogenetics/plant_tribes'
TEST_ROOT='/scratch/biotools/galaxy/test'
TOOL_DATA_TABLE=$PLANTTRIBES_TEST_ROOT/tool_data_table_conf.xml.test
WORKFLOW_ROOT='/scratch/biotools/galaxy/test/workflows'

BASE_PARAMS='--conda_auto_install --conda_debug --conda_dependency_resolution --no_cleanup'

mkdir -p $LOG_DIR

touch $LOG_FILE
echo -e "\n\n========================================" >> $LOG_FILE
date >> $LOG_FILE
echo -e "========================================\n\n" >> $LOG_FILE

for arg in "$@"; do
    case "$arg" in
        --assembly_post_processor)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "AssemblyPostProcessor\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='assembly_post_processor'
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --gene_family_classifier)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "GeneFamilyClassifier\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='gene_family_classifier'
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --gene_family_integrator)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "GeneFamilyIntegrator\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='gene_family_integrator'
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --gene_family_aligner)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "GeneFamilyAligner\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='gene_family_aligner'
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --gene_family_phylogeny_builder)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "GeneFamilyPhylogenyBuilder\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='gene_family_phylogeny_builder'
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --kaks_analysis)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "KaKsAnalysis\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='kaks_analysis'
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;

        --ks_distribution)
            echo -e "-------------------------------------\n" >> $LOG_FILE
            echo -e "KsDistribution\n" >> $LOG_FILE
            echo -e "-------------------------------------\n" >> $LOG_FILE
            TOOL_NAME='ks_distribution'
            TOOL_DIR=$PLANTTRIBES_TEST_ROOT/$TOOL_NAME
            TOOL_OUTPUTS_DIR=$LOG_DIR/$TOOL_NAME
            JOB_OUTPUT_FILES=$TOOL_OUTPUTS_DIR/job_output_files
            TEST_OUTPUT=$TOOL_OUTPUTS_DIR/$TOOL_NAME.html
            TEST_OUTPUT_JSON=$TOOL_OUTPUTS_DIR/$TOOL_NAME.json

            mkdir -p $TOOL_OUTPUTS_DIR

            echo $TOOL_NAME >> $LOG_FILE
            planemo -v lint $TOOL_DIR/$TOOL_NAME.xml >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test $BASE_PARAMS --conda_ensure_channels $CONDA_CHANNELS --conda_exec $CONDA_EXEC --conda_prefix $CONDA_PREFIX --galaxy_root $GALAXY_ROOT_PRISTINE --job_output_files $JOB_OUTPUT_FILES --test_data $TOOL_DIR/test-data --test_output $TEST_OUTPUT --test_output_json $TEST_OUTPUT_JSON $TOOL_DIR >> $LOG_FILE
            echo -e "\n\n" >> $LOG_FILE

            planemo -v test_reports $TEST_OUTPUT_JSON --test_output_text $TOOL_OUTPUTS_DIR/$TOOL_NAME.txt
            ;;
    esac
done

