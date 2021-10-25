#!/usr/bin/env python

directories = [
  '/eos/uscms/store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/SingleMuon_Overlap_4GeV/ParticleGuns/CRAB3/190125_014345/0000/',
  '/eos/uscms/store/group/l1upgrades/L1MuonTrigger/P2_10_4_0/SingleMuon_Overlap_4GeV/ParticleGuns/CRAB3/190125_014345/0001/',
]

outfile = '/tmp/jiafu/ntuple.root'

def call_cmd(cmd):
  import shlex, subprocess
  p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
  lines = p.stdout.read().split()
  return lines

def list_input_files(directories):
  all_lines = []
  for directory in directories:
    cmd = 'xrdfs root://cmseos.fnal.gov ls -u {0}'.format(directory)
    lines = call_cmd(cmd)
    lines = [line for line in lines if line.endswith('.root')]
    all_lines += lines
  return ' '.join(all_lines)

# Main
if __name__ == '__main__':
  infiles = list_input_files(directories)
  cmd = 'hadd -f {0} {1}'.format(outfile, infiles)
  lines = call_cmd(cmd)
  #print '\n'.join(lines)