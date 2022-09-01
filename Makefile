runall: test
	for example in `ls examples/*.rs`; do cargo run --example `basename $$example .rs`; done

buildall: test
	cargo build

test:
	cargo test

clean:
	$(RM) img/*.png
	$(RM) -r target

.PHONY : runall buildall clean
