LOCAL_BIN=/usr/local/bin
APP=$(LOCAL_BIN)/rvec
APP_DATA=$${HOME}/.rvec
SETTINGS=./config/config.yaml

PLATFORM:=$(shell uname)

.PHONY: test
.DEFAULT_GOAL: build

ifeq ($(PLATFORM),Windows_NT)
	$(error Windows not supported via make)
endif

install:
	mkdir -p $(APP_DATA)

	cp $(SETTINGS) $(APP_DATA)/config.yaml
	cp -r ./vendor/primer3_config $(APP_DATA) 
	cp -r ./assets/addgene/db/** $(APP_DATA) 
	cp -r ./assets/igem/db/** $(APP_DATA)
	cp -r ./assets/dnasu/db/** $(APP_DATA)
	cp ./assets/snapgene/features.tsv $(APP_DATA)
	cp ./assets/neb/enzymes.tsv $(APP_DATA)

ifeq ($(PLATFORM),Linux)
	install ./bin/linux $(APP)
	install -C ./vendor/linux/blastn $(LOCAL_BIN)
	install -C ./vendor/linux/ntthal $(LOCAL_BIN)
	install -C ./vendor/linux/primer3_core $(LOCAL_BIN)
	install -C ./vendor/linux/blastdbcmd $(LOCAL_BIN)
endif

ifeq ($(PLATFORM),Darwin)
	install ./bin/darwin $(APP)
	install -C ./vendor/darwin/blastn $(LOCAL_BIN)
	install -C ./vendor/darwin/ntthal $(LOCAL_BIN)
	install -C ./vendor/darwin/primer3_core $(LOCAL_BIN)
	install -C ./vendor/darwin/blastdbcmd $(LOCAL_BIN)
endif

build:
	rm -f ./bin/* 
	go get -d
	env GOOS=linux go build -o ./bin/linux -v
	env GOOS=darwin go build -o ./bin/darwin -v
	env GOOS=windows go build -o ./bin/rvec.exe -v

windows:
	cd scripts && makensis windows_installer.nsi

all: build install

dbs:
	cd assets && sh makeblastdbs.sh

uninstall: clean
	rm $(APP)
	rm -rf $(APP_DATA)

test:
	go test -timeout 200s ./internal/rvec
