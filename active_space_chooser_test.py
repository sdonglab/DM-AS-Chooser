from unittest import mock
import pathlib
import active_space_chooser as asc

import pytest

def test_process_opts_gdm(tmp_path):
    parser, _, _ = asc.get_parsers()
    data_dir = tmp_path / 'data'
    cmd = ['gdm-as', '--data-dir', str(data_dir), '-r', '0.5']
    opts = parser.parse_args(cmd)

    mock_error = mock.Mock(side_effect=SystemExit)
    mock_gdm_parser = mock.Mock(error=mock_error)

    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [mock.call(f'{data_dir} does not exist')]
    mock_error.reset_mock()
    
    opts = parser.parse_args(cmd)
    data_dir.mkdir()
    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [mock.call(f'did not find any multi-reference calculation files in {data_dir}')]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    touch(data_dir / 'bad-format' / 'foo.log')
    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [mock.call(f'did not find any multi-reference calculation files in {data_dir}')]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    good_log = touch(data_dir / '2-2' / 'foo.log')
    asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_count == 0
    assert opts.mr_files == [str(good_log)]
    mock_error.reset_mock()
    

def test_process_opts_edm(tmp_path):
    pass

def touch(path: pathlib.Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()

    return path