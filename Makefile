runall:
	for example in `ls examples/*.rs`; do cargo run --example `basename $$example .rs`; done

buildall:
	cargo build

clean:
	$(RM) data/*.txt img/*.png
	$(RM) -r target

.PHONY : runall buildall clean
