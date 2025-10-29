apply_cmssw_customization_steps() {
    run_cmd git cms-init

    run_cmd git cms-addpkg CalibMuon/CSCCalibration
    run_cmd git cms-addpkg CondFormats/CSCObjects
    # run_cmd git cms-addpkg Configuration/Generator
    # run_cmd git cms-addpkg Configuration/StandardSequences
    run_cmd git cms-addpkg DataFormats/CSCDigi
    run_cmd git cms-addpkg DataFormats/L1TMuon
    run_cmd git cms-addpkg EventFilter/CSCRawToDigi
    run_cmd git cms-addpkg L1Trigger/CSCTriggerPrimitives
    # run_cmd git checkout -b emtf_data_format_dev # uncomment if you want to create a new branch (make sure it does not exist yet) to avoid overwriting stuff
    # run_cmd git remote add valeria git@github.com:valeriadamante/cmssw.git # uncomment if you are not valeria and you want to work from the branch emtf_data_format
    # run_cmd git fetch valeria  # uncomment if you are not valeria and you want to work from the branch emtf_data_format
    # run_cmd git pull valeria emtf_data_format  #  uncomment if you are not valeria and you want to work from the branch emtf_data_format

    run_cmd git checkout emtf_data_format # uncomment if you are Valeria and you want to work from the branch emtf_data_format

    run_cmd mkdir GEMCSCTriggerTest
    run_cmd ln -s "$ANALYSIS_PATH/CSCSlopeFinder" GEMCSCTriggerTest/CSCSlopeFinder
}
# CalibMuon  CondFormats  DataFormats  EventFilter  GEMCSCTriggerTest  L1Trigger



CMSSW_VERSION="CMSSW_16_0_0_pre1"


action() {
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"
    local this_file_path="$this_dir/$(basename $this_file)"
    export ANALYSIS_PATH="$this_dir"
    echo "Running action with CMSSW version: $CMSSW_VERSION"

    source $ANALYSIS_PATH/GEM-CSC-trg-dev/env.sh "$this_file_path" "$CMSSW_VERSION" "$@"
    ln -fs "$ANALYSIS_PATH/GEM-CSC-trg-dev/luts" "$ANALYSIS_PATH/CSCSlopeFinder/luts"
}

action "$@"
unset -f action
unset -f apply_cmssw_customization_steps