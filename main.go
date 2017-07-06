package main

import (
  "bufio"
  "flag"
  "fmt"
  "io"
  "log"
  "os"
  "strings"
  "sync"
  "strconv"
  "regexp"
  "bytes"
)

var concurrency int = 10

type Config struct {
  inPath string
  bedPath string
  ucscChr bool
  chrIdx int
  posIdx int
}

func setup(args []string) *Config {
  config := &Config{}
  flag.StringVar(&config.inPath, "inPath", "", "The input file path (optional: default is stdin)")
  flag.StringVar(&config.bedPath, "bedPath", "", "The bed file path containing chr , pos to filter on")
  flag.BoolVar(&config.ucscChr, "ucscChr", false, "Normalize bed file to ucsc style chromosomes (chrN)")
  flag.IntVar(&config.chrIdx, "chrIdx", 0, "The chr field index")
  flag.IntVar(&config.posIdx, "posIdx", 1, "The position column index")

  // allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
  // can only run 1 such test, else, redefined flags error
  a := os.Args[1:]
  if args != nil {
    a = args
  }

  userSetPosIdx := regexp.MustCompile("-posIdx").MatchString(strings.Join(a, " "))

  flag.CommandLine.Parse(a)

  if userSetPosIdx == false {
    if regexp.MustCompile(`\.bim$`).MatchString(config.bedPath) {
      config.posIdx = 3
    }
  }

  return config
}

func init() {
  log.SetFlags(0)
}

func main() {
  config := setup(nil)

  if config.bedPath == "" {
    log.Fatal("bedPath is requried")
  }

  bedFh, err := os.Open(config.bedPath)

  if err != nil {
    log.Fatal("Couldn't open bed file")
  }
  defer bedFh.Close()

  bedReader := bufio.NewReader(bedFh)

  bedMap := readBed(config, bedReader)

  inFh := (*os.File)(nil)
  if config.inPath != "" {
    var err error

    inFh, err = os.Open(config.inPath)
    if err != nil {
      log.Fatal(err)
    }
  } else {
    inFh = os.Stdin
  }

  // make sure it gets closed
  defer inFh.Close()

  reader := bufio.NewReader(inFh)

  readFile(config, reader, bedMap, func(row string) {fmt.Println(row)})
}

func readBed(config *Config, reader *bufio.Reader) map[string]map[int]bool {
  existsMap := make(map[string]map[int]bool)
  chrIdx := config.chrIdx
  posIdx := config.posIdx 
  convertUCSC := config.ucscChr

  for {
    row, err := reader.ReadString('\n')

    if err == io.EOF {
      break
    } else if err != nil {
      log.Fatal("Couldn't read bed file", err)
    } else if row == "" {
      continue
    }

    record := strings.Split(row[:len(row) - 1], "\t")

    pos, err := strconv.Atoi(record[posIdx])

    if err != nil {
      log.Fatal("Couldn't convert pos to integer", record)
    }

    if convertUCSC && (len(record[chrIdx]) < 4 || record[chrIdx][:3] != "chr") {
      var buffer bytes.Buffer
      buffer.WriteString("chr")
      buffer.WriteString(record[chrIdx])
      record[chrIdx] = buffer.String()
    }

    _, ok := existsMap[record[chrIdx]];

    if !ok {
      existsMap[record[chrIdx]] = map[int]bool{ pos: true}
    } else {
      existsMap[record[chrIdx]][pos] = true
    }
  }

  return existsMap
}

func readFile(config *Config, reader *bufio.Reader, posMap map[string]map[int]bool,
  resultFunc func(row string)) {
  // Read buffer
  workQueue := make(chan string, 100)
  complete := make(chan bool)
  // Write buffer
  results := make(chan string, 100)
  var wg sync.WaitGroup

  endOfLineByte, numChars, headerLine, err := findEndOfLineChar(reader, "")

  if err != nil {
    log.Fatal(err)
  }

  if regexp.MustCompile("##fileformat=VCF").MatchString(headerLine) {
    // Skip VCF headers
    for {
      row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

      if err == io.EOF {
        break
      } else if err != nil {
        log.Fatal(err)
      } else if row == "" {
        // We may have not closed the pipe, but not have any more information to send
        // Wait for EOF
        continue
      }

      headerLine += row

      if row[:6] == "#CHROM" {
        break
      }
    }
  }

  fmt.Print(headerLine)

  // Read the lines into the work queue.
  go func() {
    for {
      row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

      if err == io.EOF {
        break
      } else if err != nil {
        log.Fatal(err)
      } else if row == "" {
        // We may have not closed the pipe, but not have any more information to send
        // Wait for EOF
        continue
      }

      workQueue <- row[:len(row) - numChars];
    }

    // Close the channel so everyone reading from it knows we're done.
    close(workQueue)
  }()

  wg.Add(1)
  go func() {
    defer wg.Done()
    for line := range results {
      resultFunc(line)
    }
  }()

  // Now read them all off, concurrently.
  for i := 0; i < concurrency; i++ {
    go processLine(config, posMap, config.ucscChr, workQueue, results, complete)
  }

  // Wait for everyone to finish.
  for i := 0; i < concurrency; i++ {
    <-complete
  }

  close(results)

  wg.Wait()
}

func processLine(config *Config, posMap map[string]map[int]bool, convertUCSC bool, queue chan string, results chan string, complete chan bool) {
  for line := range queue {
    record := strings.Split(line, "\t")

    pos, err := strconv.Atoi(record[1])

    if (err != nil) {
      log.Fatal("Couldn't convert position", record)
    }

    if convertUCSC && (len(record[0]) < 4 || record[0][:3] != "chr") {
      var buffer bytes.Buffer
      buffer.WriteString("chr")
      buffer.WriteString(record[0])
      record[0] = buffer.String()
    }

    if posMap[record[0]][pos] {
      results <- line
    }
  }

  // log.Println("Worker hit, missed this many times: ", hitCount, missCount)
  // Let the main process know we're done.
  complete <- true
}

func findEndOfLineChar (r *bufio.Reader, s string) (byte, int, string, error) {
  runeChar, _, err := r.ReadRune()

  if err != nil {
    return byte(0), 0, "", err
  }

  if runeChar == '\r' {
    nextByte, err := r.Peek(1)

    if err != nil {
      return byte(0), 0, "", err
    }

    if rune(nextByte[0]) == '\n' {
      //Remove the line feed
      _, _, err = r.ReadRune()

      if err != nil {
        return byte(0), 0, "", err
      }

      return nextByte[0], 2, s, nil
    }

    return byte('\r'), 1, s, nil
  }

  if runeChar == '\n' {
    return byte('\n'), 1, s, nil
  }

  s += string(runeChar)
  return findEndOfLineChar(r, s)
}