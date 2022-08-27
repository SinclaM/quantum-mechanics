BIN_NAMES = square_well_shooting_method_even square_well_shooting_method_odd

runall: buildall 
	for name in ${BIN_NAMES}; do cargo run --bin $$name; done

buildall:
	cargo build

clean:
	$(RM) data/*.txt img/*.png
	$(RM) -r target


