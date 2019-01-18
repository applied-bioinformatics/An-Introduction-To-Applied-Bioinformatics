# Biological Information <link src="2GgtNI"/>

Biological systems and computer systems are analogous in ways that may not be immediately apparent. Before we dive into using computers to study biology, let's briefly explore a relationship between the two: information processing is one of the most fundamental functions of both.

Information is a quantifiable concept, an idea that has its roots in Boolean algebra and in Claude Shannon's work on Information Theory. The most basic unit of information is the _binary digit_, or _bit_, which has two possible states. Depending on the domain, the symbols representing these two states might be `0` and `1`, `yes` and `no`, `+` and `-`, `true` and `false`, or `on` and `off`. When you answer a "yes/no" question in a conversation, or a "true/false" question on an exam, you're providing one bit of information.

To illustrate the amount of information contained in a bit, consider the game [Twenty Questions](https://en.wikipedia.org/wiki/Twenty_Questions). This game is played by two individuals, a questioner and an answerer. The answerer thinks of an object, and the questioner may ask yes/no questions to try to discover what the object is. The questioner wins the game if they discover what the object is with 20 or fewer questions; otherwise the answerer wins. When played strictly the only answers that are allowed are "yes" and "no". As the questioner asks more questions, they accumulate more information that will ideally help them determine what the object is.

Bits can also be used to represent numbers, and when they're used this way a binary mathematical system is being used. One bit can represent exactly two numbers, usually the numbers 0 and 1. For example:

 * `0` can represent the number 0
 * `1` can represent the number 1

If two bits of information are used, four numbers can be represented. For example:

 * `00` can represent the number 0
 * `01` can represent the number 1
 * `10` can represent the number 2
 * `11` can represent the number 3

If three bits of information are used, eight numbers can be represented. For example:

 * `000` can represent the number 0
 * `001` can represent the number 1
 * `010` can represent the number 2
 * `011` can represent the number 3
 * `100` can represent the number 4
 * `101` can represent the number 5
 * `110` can represent the number 6
 * `111` can represent the number 7

There is a pattern emerging here. If `n` is number of bits that you have available to represent a number, the number of numbers that you can represent is $2**n$.

Binary numbers form a _base 2_, or binary, mathematical system, because there are two possible symbols (0 and 1). We typically work in a _base 10_, or decimal, mathematical system, where there are 10 possible symbols (0, 1, 2, 3, 4, 5, 6, 7, 8, and 9). In a decimal system, we represent numbers larger than 9 by using multiple _places_: the _ones_ place, the _tens_ place, the _hundreds_ place, and so on. These are the exponents of 10: the ones place is $10**0$, the tens place is $10**1$, the hundreds place is $10**2$, and so on. When we write a decimal number with multiple places, such as 42, what we're representing is a four in the tens place plus a two in the ones place, or $4 x 10**1 + 2 x 10**0 = 42$.

In a binary system, we can also use multiple places to represent numbers larger than 1. In this case, we have our ones place ($2**0$), our _twos_ place ($2**1$), our _fours_ place ($2**2$), and so on. Thus the decimal interpretation of the number `011` is $0 x 2**2 + 1 x 2**1 + 1 x 2**0 = 3$. Try this for some of the three digit binary numbers listed above.

Computers, at their most fundamental level, work on the basis of binary numbers, which is why you have probably that computers represent data as zeros and ones. Thankfully, there are multiple layers of translation (like the one we just used to convert the binary number `011` to a decimal value) applied to those zeros and ones before a computer's data is presented to us humans!

_Information_ is technically defined as a sequence of symbols that can be interpreted as a message. To put these terms in the context of our examples above, our _message_ is a decimal number, our _symbols_ are 0 and 1, and the sequence is the ordered collection of symbols, such as `011`. The number of places (let's call that $p$), and the number of symbols (let's call that $n_symbols$) define the number of different messages ($n_messages$) that can be encoded as: $n_messages = n_symbols**p$.
