all:
	@echo 'RNAstructure repository testing started...'
	@echo
	make AllSub bifold bipartition ct2dot dot2ct DotPlots draw DuplexFold dynalign efn2 Fold MaxExpect NAPSS oligoscreen oligowalk PARTS partition ProbablePair refold RemovePseudoknots stochastic cleanup
	@echo 'RNAstructure repository testing finished.'

AllSub: tests/AllSub
tests/AllSub:
	@echo 'AllSub testing started...'
	@./NixScript AllSub
	@echo 'AllSub testing finished.'
	@echo

bifold: tests/bifold
tests/bifold:
	@echo 'bifold testing started...'
	@./NixScript bifold
	@echo 'bifold testing finished.'
	@echo

bipartition: tests/bipartition
tests/bipartition:
	@echo 'bipartition testing started...'
	@./NixScript bipartition
	@echo 'bipartition testing finished.'
	@echo

ct2dot: tests/ct2dot
tests/ct2dot:
	@echo 'ct2dot testing started...'
	@./NixScript ct2dot
	@echo 'ct2dot testing finished.'
	@echo

dot2ct: tests/dot2ct
tests/dot2ct:
	@echo 'dot2ct testing started...'
	@./NixScript dot2ct
	@echo 'dot2ct testing finished.'
	@echo


DotPlots: tests/DotPlots
tests/DotPlots:
	@echo 'DotPlots testing started...'
	@./NixScript DotPlots
	@echo 'DotPlots testing finished.'
	@echo

draw: tests/draw
tests/draw:
	@echo 'draw testing started...'
	@./NixScript draw
	@echo 'draw testing finished.'
	@echo

DuplexFold: tests/DuplexFold
tests/DuplexFold:
	@echo 'DuplexFold testing started...'
	@./NixScript DuplexFold
	@echo 'DuplexFold testing finished.'
	@echo

dynalign: tests/dynalign
tests/dynalign:
	@echo 'dynalign testing started...'
	@./NixScript dynalign
	@echo 'dynalign testing finished.'
	@echo

efn2: tests/efn2
tests/efn2:
	@echo 'efn2 testing started...'
	@./NixScript efn2
	@echo 'efn2 testing finished.'
	@echo

Fold: tests/Fold
tests/Fold:
	@echo 'Fold testing started...'
	@./NixScript Fold
	@echo 'Fold testing finished.'
	@echo

MaxExpect: tests/MaxExpect
tests/MaxExpect:
	@echo 'MaxExpect testing started...'
	@./NixScript MaxExpect
	@echo 'MaxExpect testing finished.'
	@echo

NAPSS: tests/NAPSS
tests/NAPSS:
	@echo 'NAPSS testing started...'
	@./NixScript NAPSS
	@echo 'NAPSS testing finished.'
	@echo

oligoscreen: tests/oligoscreen
tests/oligoscreen:
	@echo 'oligoscreen testing started...'
	@./NixScript oligoscreen
	@echo 'oligoscreen testing finished.'
	@echo

oligowalk: tests/oligowalk
tests/oligowalk:
	@echo 'oligowalk testing started...'
	@./NixScript oligowalk
	@echo 'oligowalk testing finished.'
	@echo

PARTS: tests/PARTS
tests/PARTS:
	@echo 'PARTS testing started...'
	@./NixScript PARTS
	@echo 'PARTS testing finished.'
	@echo

partition: tests/partition
tests/partition:
	@echo 'partition testing started...'
	@./NixScript pfunction
	@echo 'partition testing finished.'
	@echo

ProbablePair: tests/ProbablePair
tests/ProbablePair:
	@echo 'ProbablePair testing started...'
	@./NixScript ProbablePair
	@echo 'ProbablePair testing finished.'
	@echo

refold: tests/refold
tests/refold:
	@echo 'refold testing started...'
	@./NixScript refold
	@echo 'refold testing finished.'
	@echo

RemovePseudoknots: tests/RemovePseudoknots
tests/RemovePseudoknots:
	@echo 'RemovePseudoknots testing started...'
	@./NixScript RemovePseudoknots
	@echo 'RemovePseudoknots testing finished.'
	@echo

stochastic: tests/stochastic
tests/stochastic:
	@echo 'stochastic testing started...'
	@./NixScript stochastic
	@echo 'stochastic testing finished.'
	@echo

cleanup: tests/cleanup
tests/cleanup:
	@./NixScript cleanup

runTests: all
testing: all
tests: all
