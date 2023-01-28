from unittest import mock
import pathlib
import active_space_chooser as asc

import pytest

class TestGDMSelector:
    @mock.patch.object(asc.GDMSelector, 'get_ground_state_dipole')
    def test_select(self, mock_get_dipole):
        mock_get_dipole.side_effect = [0.1, 0.3, 0.2]
        mr_calcs = [
            asc.MultiRefCalc(2, 2, 'foo'),
            asc.MultiRefCalc(4, 4, 'bar'),
            asc.MultiRefCalc(6, 6, 'baz'),
        ]
        selector = asc.GDMSelector(mr_calcs=mr_calcs, ref_dipole=0.12)
        assert selector.select() == mr_calcs[0], "gdm selects closest to ref"

        mock_get_dipole.side_effect = [0.1, 0.14, 0.2]
        mr_calcs = [
            asc.MultiRefCalc(2, 2, 'foo'),
            asc.MultiRefCalc(4, 4, 'bar'),
            asc.MultiRefCalc(6, 6, 'baz'),
        ]
        selector = asc.GDMSelector(mr_calcs=mr_calcs, ref_dipole=0.12)
        assert selector.select() == mr_calcs[0], "gdm select first active space found in case of ties"

class TestEDMSelector:
    def test_select(self):
        pass


def test_process_opts_gdm(tmp_path):
    parser, _, _ = asc.get_parsers()
    data_dir = tmp_path / 'data'
    cmd = ['gdm-as', '--data-dir', str(data_dir), '-r', '0.5']
    opts = parser.parse_args(cmd)

    mock_error = mock.Mock(side_effect=SystemExit)
    mock_gdm_parser = mock.Mock(error=mock_error)

    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [
        mock.call(f'{data_dir} does not exist')
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    data_dir.mkdir()
    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [
        mock.call(
            f'did not find any multi-reference calculation files in {data_dir}'
        )
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    touch(data_dir / 'bad-format' / 'foo.log')
    with pytest.raises(SystemExit):
        asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_args_list == [
        mock.call(
            f'did not find any multi-reference calculation files in {data_dir}'
        )
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    foo_log = touch(data_dir / '2-2' / 'foo.log')
    asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_count == 0
    assert opts.mr_files == [str(foo_log)]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    bar_log = touch(data_dir / '3-3' / 'bar.log')
    asc.process_opts(mock_gdm_parser, None, opts)
    assert mock_error.call_count == 0
    assert opts.mr_files == [str(foo_log), str(bar_log)]
    mock_error.reset_mock()


def test_process_opts_edm(tmp_path):
    parser, _, _ = asc.get_parsers()
    data_dir = tmp_path / 'data'
    cmd = ['edm-as', '--data-dir', str(data_dir)]
    opts = parser.parse_args(cmd)

    mock_error = mock.Mock(side_effect=SystemExit)
    mock_edm_parser = mock.Mock(error=mock_error)

    with pytest.raises(SystemExit):
        asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_args_list == [
        mock.call(f'{data_dir} does not exist')
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    data_dir.mkdir()
    with pytest.raises(SystemExit):
        asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_args_list == [
        mock.call(
            f'did not find any multi-reference calculation files in {data_dir}'
        )
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    touch(data_dir / 'bad-format' / 'foo.log')
    with pytest.raises(SystemExit):
        asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_args_list == [
        mock.call(
            f'did not find any multi-reference calculation files in {data_dir}'
        )
    ]
    mock_error.reset_mock()

    opts = parser.parse_args(cmd)
    foo_log = touch(data_dir / '2-2' / 'foo.log')
    with pytest.raises(SystemExit):
        asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_args_list == [
        mock.call(f'did not find td-dft calculation file in {data_dir}')
    ]
    mock_error.reset_mock()

    tddft_log = data_dir / 'test.log'
    cmd = ['edm-as', '--data-dir', str(data_dir), '-r', str(tddft_log)]
    opts = parser.parse_args(cmd)
    with pytest.raises(SystemExit):
        asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_args_list == [
        mock.call(f'{tddft_log} does not exist')
    ]
    mock_error.reset_mock()

    touch(tddft_log)

    opts = parser.parse_args(cmd)
    asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_count == 0
    assert opts.mr_files == [str(foo_log)]
    assert opts.ref_tddft == str(tddft_log)
    mock_error.reset_mock()

    cmd = ['edm-as', '--data-dir', str(data_dir)]
    asc.process_opts(None, mock_edm_parser, opts)
    assert mock_error.call_count == 0
    assert opts.mr_files == [str(foo_log)]
    assert opts.ref_tddft == str(tddft_log)
    mock_error.reset_mock()


def touch(path: pathlib.Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()

    return path
