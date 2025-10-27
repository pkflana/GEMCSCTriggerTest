apply_cmssw_customization_steps() {
    run_cmd git cms-addpkg DataFormats/LCTDebug
    run_cmd git cms-addpkg DataFormats/CSCDigi
    run_cmd git cms-addpkg DataFormats/L1TMuon
    run_cmd git cms-addpkg Configuration/Generator
    run_cmd git cms-addpkg Configuration/StandardSequences
    run_cmd git cms-addpkg CondFormats/CSCObjects
    run_cmd git cms-addpkg CalibMuon/CSCCalibration
    run_cmd git cms-addpkg EventFilter/CSCRawToDigi
    run_cmd mkdir GEMCSCTriggerTest
    run_cmd ln -s "$ANALYSIS_PATH/CSCSlopeFinder" GEMCSCTriggerTest/CSCSlopeFinder
}


CMSSW_VERSION="CMSSW_14_2_0_pre1"

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