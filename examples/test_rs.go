package main

import (
	"os"

	"fmt"

	"github.com/klauspost/reedsolomon"
)

func main() {
	dataShards := 12
	parShards := 6
	enc, err := reedsolomon.New(dataShards, parShards)
	checkErr(err)

	text := "hello world hello world "
	data := make([]byte, len(text))
	copy(data, []byte(text))
	fmt.Println("data len:", len(data))

	shards, err := enc.Split(data)
	checkErr(err)
	fmt.Println("split shards:", shards)

	err = enc.Encode(shards)
	checkErr(err)
	fmt.Println("encode shards:", shards)

	dec, err := reedsolomon.New(dataShards, parShards)
	checkErr(err)
	//fmt.Println(len(shards), len(shards[0]))

	shards[1] = nil
	shards[3] = nil
	shards[4] = nil
	//fmt.Println(shards)

	err = dec.Reconstruct(shards)
	checkErr(err)
	fmt.Println("decode:", shards)
}

func checkErr(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err.Error())
		os.Exit(2)
	}
}
