from __future__ import annotations

from typing import Optional


class Trinket:
    """A useless object that stores made-up data for no scientific purpose.

    Attributes:
        name: A whimsical name for the trinket.
        shininess: How shiny the trinket is, on a scale from 0 to 1.
        num_sides: The number of sides the trinket has.
        secret_message: An optional cryptic message hidden inside.
    """

    def __init__(
        self,
        name: str,
        shininess: float = 0.5,
        num_sides: int = 7,
        secret_message: Optional[str] = None,
    ) -> None:
        """Initializes a Trinket with made-up properties.

        Args:
            name: A whimsical name for the trinket.
            shininess: How shiny the trinket is, from 0 (dull) to 1 (blinding).
            num_sides: The number of sides the trinket has. Seven is the default
                because seven is a lucky number.
            secret_message: An optional cryptic message hidden inside.
        """
        self.name = name
        self.shininess = shininess
        self.num_sides = num_sides
        self.secret_message = secret_message

    def describe(self) -> str:
        """Returns a human-readable description of the trinket.

        Returns:
            A sentence describing the trinket's properties.
        """
        description = (
            f"'{self.name}' is a {self.num_sides}-sided trinket "
            f"with a shininess of {self.shininess:.2f}."
        )
        if self.secret_message is not None:
            description += f" It whispers: '{self.secret_message}'."
        return description

    def is_blinding(self) -> bool:
        """Returns True if the trinket's shininess exceeds safe levels.

        Returns:
            True if shininess is above 0.9, False otherwise.
        """
        return self.shininess > 0.9
