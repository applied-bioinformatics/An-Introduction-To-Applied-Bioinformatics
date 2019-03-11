# Biological Information <link src="2GgtNI"/>

Biological systems and computer systems are analogous in ways that may not be immediately apparent. Before we dive into using computers to study biology, let's briefly explore a relationship between the two: information processing is one of the most fundamental functions of both.

## Central Dogma of Molecular Biology

The Central Dogma of Molecular Biology describes information flow in biological systems. It begins with DNA, a relatively long-lived information storage molecule, from which information typically flows in two directions: into new DNA molecules, during the process of replication, or into messenger RNA (mRNA), during the processing of transcription. mRNA is a relatively short-lived molecule that transfers information that is used to synthesize protein molecules by the ribosome. Proteins are often thought of as the building blocks of life. They serve a variety of purposes, ranging from molecular machines such as transmembrane ion transporters, to structural molecules like myosin, a major component of muscle fibers. There are some uncommon circumstances where information flows differently, for example some viruses can reverse transcribe RNA to DNA, but proteins seem to be a terminal output of this information flow. Once a protein has been created, we are aware of no process that can work backwards to re-create the RNA or DNA that encoded it.

We'll revisit these ideas at the end of this chapter, but first let's establish some concepts that will help us to understand and even quantify information. These ideas have their roots in Boolean algebra and Information Theory. Bear with me while I introduce some concepts that may be new to you, and may initially seem unrelated.

TODO: creative commons central dogma figure.

## Binary and decimal numerical systems

Humans most frequently use a _base 10_ or decimal numerical system for representing numbers. _Base 10_ means than there are ten available digits including zero. These are the digits 0, 1, 2, 3, 4, 5, 6, 7, 8, and 9. We represent numbers larger than 9 using multiple places: the _ones_ place, the _tens_ place, the _hundreds_ place, and so on. These are the exponents of 10: the ones place is $10**0$, the tens place is $10**1$, the hundreds place is $10**2$, and so on. When we write a decimal number with multiple places, such as 42, what we're representing is a four in the tens place plus a two in the ones place, or $4 x 10**1 + 2 x 10**0 = 42$.

You've probably heard that computers use a _base 2_ or binary numerical system to represent numbers. The _base_ again describes the number of available digits, so in a base 2 or binary system, there are two digits, 0 and 1. These are defined as the binary digits. As in the decimal system, numbers larger than 1 are represented using multiple places. The places in a binary number are again based on exponents, but this time they are the exponents of 2. Instead of a ones place, a tens place, and a hundreds place, the first three places in a binary number are the ones place ($2**0$), the _twos_ place ($2**1$), and the _fours_ place ($2**2$). Thus the interpretation of the binary number `011` is $0 x 2**2 + 1 x 2**1 + 1 x 2**0 = 3$.

When working with numbers that may be other than base 10, by convention numbers would be written as $(n)_b$, where $n$ is the number, and $b$ is the base of the number. For example, $(11)_10$ represents the decimal number 11, because the base is 10. $(11)_2$ represents the decimal number 3: because the base is 2, we know that this is a binary number.

Here are some binary numbers and formulas for translating them to their decimal equivalents.

 * $(0)_2$ is the decimal number 0 ($0 x 2**0$)
 * $(1)_2$ is the decimal number 1 ($1 x 2**0$)
 * $(01)_2$ is also the decimal number 1 ($0 x 2**1 + 1 x 2**0$)
 * $(11)_2$ is the decimal number 3 ($1 x 2**1 + 1 x 2**0$)
 * $(110)_2$ is the decimal number 6 ($1 x 2**2 + 1 x 2**1 + 0 x 2**0$)
 * $(111)_2$ is the decimal number 7 ($1 x 2**2 + 1 x 2**1 + 1 x 2**0$)

A single **bi**nary digi**t** (a zero or one) is referred to as a _bit_, and bits can be used to encode a lot more than just numbers.

## Encoding messages in bits

The messages encoded by bits can be nearly anything, provided that the sender of the message and the recipient of the message have agreed on a coding scheme which describes how a message can be encoded in bits or decoded from bits. The number of messages that can be sent using bits is a simple function of the number of places that are used. For example, if the only messages I want to transmit to you are "yes" and "no", I could achieve that by transmitting a single bit of information. You and I could agree that "yes" will be represented by the bit 1, and "no" will be represented by the bit 0.

Internally, computers send and receive messages that are encoded using electrical currents. To reduce errors in message transmission, the electrical currents are interpreted only as being off or on, such that a message may be transmitted as off-on-on. These two states, off and on, are often interpreted by computers as binary numbers, where zero is synonymous with off (no current) and one is synonymous with on (current). So our message of off-on-on could be read as the binary number 011, or the decimal number 3.

To illustrate a useful system that operates on the transmission of one bit of information, I'll describe a photosensor for an outdoor spotlight. In this example the photosensor is the sender of the message and the spotlight is the receiver of the message. The transmission of a zero from the photosensor to the spotlight could mean that it is currently light outside, and the transmission of a one could mean that it is currently dark outside. (The meanings of zero and one could be reversed: all that matters is that the sender and the receiver know what each value means.) The photosensor can monitor the available light, and send a message to the spotlight once per minute. If it is currently light outside, the photosensor will send a zero to the spotlight and the spotlight will turn off or remain off. If it is currently dark outside, the photosensor will send a one to the spotlight, and the spotlight will turn on or remain on. The photosensor is functioning as an on/off switch for the spotlight, transmitting one bit of information every minute.

There are couple of important things to consider in this example. First, the meaning of "currently light outside" and "currently dark outside" are embodied in the photosensor. It must make a decision on whether it is light or dark on it's own, because it is only transmitting one bit of information (zero equals light and one equals dark). The message it sends isn't complex enough to describe how light or dark it currently is outside -it's effectively only flipping a switch on or off.

To enable the transmission of more complex messages more bits can be used. One bit allows us to transmit two messages: 0 and 1, which in our photosensor example are interpreted as _off_ and _on_, respectively. If our message is based on two bits we can transmit four messages, 00, 01, 10, or 11. A real-world example of this could be a light switch with four states: off, low brightness, medium brightness, and high brightness. If our message is based on three bits, we can transmit eight messages, 000, 001, 010, 011, 100, 101, 110, or 111. There is a pattern emerging here. If `n` is number of bits that you have available to send a message, the number of distinct messages that you can send is $2**n$. To generalize this formula further, if the number of available digits in the system is `b`, and the number of places you can use in your message is `n`, then the number of messages that can be sent is $b**n$.

In computer systems, the bit is the most fundamental unit of information. The next largest unit is the byte, which is composed of eight bits. How many messages can be encoded in one byte?

## DNA messages are encoded in a base 4 systems

The building blocks of DNA are four chemical compounds called adenine, cytosine, guanine, and thymine. We often represent these compounds with the abbreviations A, C, G and T, respectively. Like computer systems, information is represented based on distinct states, but in biological systems there are four states rather than two. Each position or place in a DNA sequence can contain one of these compounds, and the linear order of the compounds can encode a message.

One type of message that can be encoded in DNA is the sequence of amino acids in a protein. When first translated, proteins are composed of simpler units, the amino acids, and most organisms use 20 different amino acids to build proteins. Because there are four DNA bases (A, C, G, and T) and twenty amino acids, we need more than one base to transmit the message of what amino acid comes next in a protein from DNA to the ribosome. How many DNA bases we need depends on how many messages we want to be able to send, which in this case is 20 (for the twenty amino acids). So, how many DNA bases are needed to encode the 20 canonical amino acids?

As mentioned above, we can determine the number of messages we can send in a given numerical system with a given number of places using the formula $b**n$. For messages encoded in DNA, `b` is four, so with one place (or one DNA base) we can send four messages. Since four is less than twenty, we'll need longer messages to encode the twenty amino acids. If our message were composed of two bases, we could send $4**2=16$ messages - that's still less than twenty, so we'll need more bases. If our message were composed of three bases, we could send $4**3=64** messages. This is more than twenty, which means that we can encode all of the amino acids (with some messages to spare) in three bases. It's important to note that the number of places we can use must be a whole number - "2.5 bases of DNA" is not a meaningful quantity.

Amino acids are in fact encoded by three nucleotide bases, and the three base messages are referred to as _codons_. The mapping of codons to amino acids is referred to as the _genetic code_. Each codon represents exactly one amino acid, with the exception of some, the _stop codons_, which indicate the end of a message. Because there are 64 codons but only twenty-one messages that need to be transmitted (the twenty amino acids and the "stop" signal), some amino acids and the stop signal are represented by more than one codon. This is referred to as the redundancy of the genetic code.

<figure>
    <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/getting-started/images/genetic-code.png">
    <figcaption>Figure 2: The vertebrate RNA genetic code. The corresponding DNA genetic code is identical, except that Us are replaced with Ts. (Figure attribution: NIH [Public domain], <a href="https://commons.wikimedia.org/wiki/File:06_chart_pu3.png">via Wikimedia Commons</a>.)</figcaption>
</figure>
<p>

TODO: simple skbio translation example

## Information theory

[Introduce entropy and mutual information, tie back to information loss in RNA to protein]

Information is a quantifiable concept, an idea that has its roots in Boolean algebra and in Claude Shannon's work on Information Theory. The most basic unit of information is the _binary digit_, or _bit_, which has two possible states. Depending on the domain, the symbols representing these two states might be `0` and `1`, `yes` and `no`, `+` and `-`, `true` and `false`, or `on` and `off`. When you answer a "yes/no" question in a conversation, or a "true/false" question on an exam, you're providing one bit of information.

_Information_ is technically defined as a sequence of symbols that can be interpreted as a message. To put these terms in the context of our examples above, our _message_ is a decimal number, our _symbols_ are 0 and 1, and the sequence is the ordered collection of symbols, such as `011`. The number of places (let's call that $p$), and the number of symbols (let's call that $n_symbols$) define the number of different messages ($n_messages$) that can be encoded as: $n_messages = n_symbols**p$.
